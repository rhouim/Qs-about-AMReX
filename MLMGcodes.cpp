    // ......
    // Only consider the single level
    const int c_lev = 0;
    const int nlevel = 1;
       
    // Initialize mg_geom, mg_grids, mg_dmap
    // AMREX_SPACEDIM = 2
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        
	if (Geometry::isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            if (phys_bc->lo(idim) == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc->lo(idim) == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                // In my situation, all bc types are Neumann
		mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (phys_bc->hi(idim) == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc->hi(idim) == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                // In my situation, all bc types are Neumann
               mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    Vector<Geometry> mg_geom(nlevel);
    for (int lev = 0; lev < nlevel; lev++) {
        mg_geom[lev] = parent->Geom(lev+c_lev); 
    }  

    Vector<BoxArray> mg_grids(nlevel);
    for (int lev = 0; lev < nlevel; lev++) {
        mg_grids[lev] = parent->boxArray(lev+c_lev);
    }

    Vector<DistributionMapping> mg_dmap(nlevel);
    for (int lev=0; lev < nlevel; lev++ ) {
        mg_dmap[lev] = LevelData[lev+c_lev]->get_new_data(State_Type).DistributionMap();
    }

    // Build LPInfo and mlabec
    LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    MLABecLaplacian mlabec(mg_geom, mg_grids, mg_dmap, info);
    mlabec.setMaxOrder(2);

    // SetDomainBC 
    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);

    // Set Solu and level BC
    // phi_rebase represents the solution
    Vector<MultiFab*> phi_rebase(phi.begin()+c_lev, phi.begin()+c_lev+nlevel);
    for (int ilev = 0; ilev < nlevel; ++ilev)
    {
	    mlabec.setLevelBC(ilev, phi_rebase[ilev]);
    }

    // Assign memory for acoef, bcoef, rhs
       acoef.resize(nlevel);
       bcoef.resize(nlevel);
       rhs.resize(nlevel);
   
    // define acoef, bcoef, rhs 
    // initialize acoef as 1, rhs as 0
    for (int ilev = 0; ilev < nlevel; ++ilev)
    {
        acoef[ilev].define(mg_grids[ilev], mg_dmap[ilev], 1, 0);
        // 1 ghost cell for bcoef
	      bcoef[ilev].define(mg_grids[ilev], mg_dmap[ilev], 1, 1);
	      rhs[ilev].define(mg_grids[ilev], mg_dmap[ilev], 1, 0);

	      acoef[ilev].setVal(1,0,1);
	      rhs[ilev].setVal(0,0,1);
    }
    
    // Set scalars
    mlabec.setScalars(0.0,-1.0);

    // Copy the sig to the bcoef
    // *sig[lev] is a known mf
    for (int ilev = 0; ilev < nlevel; ++ilev) {
       MultiFab::Copy(bcoef[c_lev+ilev],*sig[c_lev+ilev],0,0,1,1);
    }

    //Set a,b coefficients of the object mlabec
        for (int ilev = 0; ilev < nlevel; ++ilev)
        {            
            mlabec.setACoeffs(ilev, acoef[ilev]);

            std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face({AMREX_D_DECL(&face_bcoef[0],
                                                            &face_bcoef[1],
                                                            &face_bcoef[2])},
                                               bcoef[ilev], mg_geom[ilev]);

            mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
        }

    //Build MLMG
    MLMG mlmg(mlabec);
    mlmg.setMaxIter(100);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(2);
    mlmg.setCGVerbose(0);
    
    // Step 1
    // Initialize the edge centered velocity mf velec_in
    Vector< std::array<MultiFab,AMREX_SPACEDIM> > velec_in(nlevel);
    // Initialize its ba and dm
       for (int ilev = 0; ilev < nlevel; ++ilev)
        { 
	      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
             const BoxArray& fluxba = amrex::convert(mg_grids[ilev], IntVect::TheDimensionVector(idim));
             (velec_in[ilev][idim]).define(fluxba, mg_dmap[ilev], 1, 0);
             (velec_in[ilev][idim]).setVal(0.0);
            }
        }

    // Step 2, calculate velec and rhs
    Vector<BCRec> bc(AMREX_SPACEDIM); // velocities bc as inputs
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
    	bc[idim]    = NavierStokesBase::get_desc_lst()[State_Type].getBC(idim);
    }

    Vector<MultiFab*> vel_in{vel.begin()+c_lev, vel.begin()+c_lev+nlevel};

    // Input: 
    // vel_in: cell ccentered velocities, known mf with 1 ghost cell
    // Outputs:
    // 1. edge centered velocity mf velec_in
    // 2. rhs: right handside for MLMG solver
    mlabec.compRHS2d(amrex::GetVecOfPtrs(rhs), velec_in, vel_in, bc);

    // check rhs
    // Here is how I get the rhs for CPU 0 and CPU 1 
        for (int ilev = 0; ilev < nlevel; ++ilev)
        { 
     		  for (MFIter vmfi(rhs[ilev]); vmfi.isValid(); ++vmfi)
    		  {
          		amrex::Print() << " Multi_method: Step 2, check rhs " << "\n";
         	  	amrex::Print() << "ilev " << ilev << "\n";
              amrex::Print() << "rhs, CPU0 " << "\n";
			        amrex::Print() << (rhs[ilev])[vmfi] << "\n";
			        amrex::Print() << "\n";
              amrex::Print(1) << "rhs, CPU1 " << "\n";
			        amrex::Print(1) << (rhs[ilev])[vmfi] << "\n";
			        amrex::Print(1) << "\n";
    		  }
        }

    // Step 3
    // solve
    // Here, if I use 2 CPUs, mlmg can not converge
    Real mlmg_err_in = mlmg.solve(phi_rebase, amrex::GetVecOfConstPtrs(rhs), 1.0e-12, 1.0e-16);
    

