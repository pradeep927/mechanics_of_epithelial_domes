
// C++ headers
#include <iostream>
#include <mpi.h>

// Trilinos headers
#include <Teuchos_RCP.hpp>
#include <math.h>

// hiperlife headers
#include "hl_TypeDefs.h"
#include "hl_Geometry.h"
#include "hl_StructMeshGenerator.h"
#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"
#include "hl_ConsistencyCheck.h"
#include <hl_LoadBalance.h>
#include <hl_LinearSolver_Direct_Amesos2.h>
#include <hl_NonlinearSolver_NewtonRaphson.h>
#include <hl_MeshLoader.h>
#include <hl_ConfigFile.h>
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <hl_LocMongeParam.h>



#include "hl_Parser.h"
#include "hl_ParamStructure.h"
#include "hl_Core.h"

// Header to auxiliary functions
#include "AuxViscousInextMembrane.h"



int main (int argc, char *argv[])
{
    using namespace std;
    using namespace hiperlife;
    using Teuchos::rcp;
    using Teuchos::RCP;
    using namespace hiperlife::Tensor;


    // **************************************************************//
    // *****                 INITIALIZATION                     *****//
    // **************************************************************//
    hiperlife::Init(argc, argv);
    const int myRank = hiperlife::MyRank();

    // Model parameters and solver options are read from a configuration file
    if (argc < 2)
    {
        cout << "Configuration file not provided!" << endl;
        MPI_Finalize();
        return 1;
    }
    // const char *config_filename = argv[1];
    ConfigFile config(argv[1]);


    //Cortex model parameters
    int aSteps{100},control_fric{10000};
    double  R{},thick{},Ainc{},aap{1.0},kappa{0.0},bbn{1.0},width{40.0},bound{0.0},bc_const{1.0},ang_old{0.0},mesh_refine{0.0},mom_steps{100.0},kspr{0.0},sp_gap{0.008},a11{0.51},a12{0.49},a33{0.98},width_kap{40.0},kap1{1.0};
    {
        config.readInto(R, "R");
        config.readInto(thick, "thick");
        config.readInto(aap, "aap");
        config.readInto(kappa, "kappa");
        config.readInto(aSteps, "aSteps");
        config.readInto(bbn,"bbn");
        config.readInto(width,"width");
        config.readInto(bound,"bound");
        config.readInto(bc_const,"bc_const");
        config.readInto(control_fric,"control_fric");
        config.readInto(ang_old,"ang_old");
        config.readInto(mesh_refine,"mesh_refine");
        config.readInto(mom_steps,"mom_steps");
        config.readInto(kspr,"kspr");
        config.readInto(sp_gap,"sp_gap");

        config.readInto(a11,"a11");
        config.readInto(a12,"a12");
        config.readInto(a33,"a33");
        config.readInto(width_kap,"width_kap");
        config.readInto(kap1, "kap1");

    }

    double mu{10.0}, fric{0.01},  young{25.0}, poisson{0.25},force{0.0},Kconf{0.0}, gamma{0.0},f0{0.9},g_ratio{0.8},v_target{0.667},vstep{0.0}, v_r_step{0.0},coeff{0.0},ang{0.0},tens_factor{1.0},out_choice{0.0};
    {
        config.readInto(fric, "fric");
        config.readInto(young, "young");
        config.readInto(poisson, "poisson");
        config.readInto(force, "force");
        config.readInto(Kconf, "Kconf");
        config.readInto(gamma, "gamma");
        config.readInto(ang, "ang");
        config.readInto(tens_factor,"tens_factor");
        config.readInto(f0, "f0");
        config.readInto(g_ratio, "g_ratio");
        config.readInto(v_target, "v_target");
        config.readInto(vstep, "vstep");
        config.readInto(v_r_step, "v_r_step");
        config.readInto(coeff, "coeff");
        config.readInto(out_choice, "out_choice");

    }

    v_target=v_target*(R);
    //Simulation time settings
    int nSave{},spring{0},hl_print{1000};
    double deltat{}, deltatMax{}, totalTime{}, tSimu{}, stepFactor{},factor{1.0},fric_fact{1},height_in{0.05},fric2{0.015}, v_rn_step{1000.0}, vol_inc{1.0},control_dt{0.0001},v_cont{0.98},out_fact{1.2},thin_shell{0.0};
    int restart{}, totalSteps{}, nPrint{1},ConsCheck{},nstep{50} ,vnstep{500},tens_start{2},choice{100}, control_steps{400}, gap{2000}, test_jac{0};
    {
        config.readInto(deltat, "deltat");
        config.readInto(deltatMax, "deltatMax");
        config.readInto(totalTime, "totalTime");
        config.readInto(totalSteps, "totalSteps");
        config.readInto(tSimu, "tSimu");
        config.readInto(stepFactor, "stepFactor");
        config.readInto(restart, "restart");
        config.readInto(nPrint, "nPrint");
        config.readInto(ConsCheck, "ConsCheck");
        config.readInto(nPrint, "nPrint");
        config.readInto(nstep, "nstep");
        config.readInto(vnstep, "vnstep");
        config.readInto(factor, "factor");
        config.readInto(tens_start, "tens_start");
        config.readInto(fric_fact, "fric_fact");
        config.readInto(height_in, "height_in");
        config.readInto(v_rn_step, "v_rn_step");
        config.readInto(v_cont, "v_cont");

        config.readInto(out_fact, "out_fact");
        config.readInto(hl_print, "hl_print");

        config.readInto(thin_shell, "thin_shell");



        config.readInto(choice, "choice");
        config.readInto(fric2, "fric2");
        config.readInto(spring, "spring");

        config.readInto(gap, "gap");
        config.readInto(test_jac, "test_jac");

        config.readInto(control_steps, "control_steps");
        config.readInto(control_dt, "control_dt");

    }



    //Numerical parameters
    double  gPts{};
    int MAXITER_AZ{}, MAXITER_NR{};
    double SOLTOL_NR{}, RESTOL_NR{}, SOLTOL_AZ{};
    string prefixMesh;
    {
        config.readInto(gPts, "gPts");
        config.readInto(MAXITER_NR, "MAXITER_NR");
        config.readInto(SOLTOL_NR, "SOLTOL_NR");
        config.readInto(RESTOL_NR, "RESTOL_NR");
        config.readInto(MAXITER_AZ, "MAXITER_AZ");
        config.readInto(SOLTOL_AZ, "SOLTOL_AZ");
        config.readInto(prefixMesh, "prefixMesh");
    }


    //Set structure

    SmartPtr<ParamStructure> paramStr = CreateParamStructure<MembParams>();
    paramStr->setRealParameter(MembParams::deltat, deltat);

    paramStr->setRealParameter(MembParams::fric, fric);

    paramStr->setRealParameter(MembParams::young, young);

    paramStr->setRealParameter(MembParams::poisson, poisson);

    paramStr->setRealParameter(MembParams::thick, thick);

    paramStr->setRealParameter(MembParams::force, force);
    paramStr->setRealParameter(MembParams::gamma, gamma);
    paramStr->setIntParameter(MembParams::tens_start, tens_start);

    paramStr->setRealParameter(MembParams::fric_fact, fric_fact);
    paramStr->setRealParameter(MembParams::fric2, fric2);
    paramStr->setRealParameter(MembParams::aap, aap);
    paramStr->setRealParameter(MembParams::R, R);

    paramStr->setRealParameter(MembParams::kappa, kappa);
    paramStr->setRealParameter(MembParams::bbn, bbn);

    paramStr->setRealParameter(MembParams::width, width);
    paramStr->setRealParameter(MembParams::tens_factor, tens_factor);
    paramStr->setIntParameter(MembParams::control_fric, control_fric);

    paramStr->setRealParameter(MembParams::height_in,height_in);
    paramStr->setRealParameter(MembParams::Kconf, Kconf);


    paramStr->setRealParameter(MembParams::f0,f0);


    paramStr->setRealParameter(MembParams::g_ratio, g_ratio);
    paramStr->setRealParameter(MembParams::bc_const,bc_const);
    paramStr->setRealParameter(MembParams::ang, ang);
    paramStr->setRealParameter(MembParams::ang_old,ang_old);
    paramStr->setRealParameter(MembParams::kspr, kspr);
    paramStr->setIntParameter(MembParams::spring,spring);
    paramStr->setRealParameter(MembParams::sp_gap, sp_gap);
    paramStr->setRealParameter(MembParams::a11,a11);
    paramStr->setRealParameter(MembParams::a12,a12);
    paramStr->setRealParameter(MembParams::a33,a33);
    paramStr->setRealParameter(MembParams::kap1,kap1);
    paramStr->setRealParameter(MembParams::width_kap,width_kap);

    paramStr->setIntParameter(MembParams::fric_start,vnstep+gap);




        //Output
        if (myRank == 0)
        {
            cout << endl << "Parameters: " << endl;
            cout << "Radius:   " << R << endl;
            cout << "thickness " <<height_in << endl;
            cout << "spring choice: 0-no, 1: yes: " <<spring << endl;

            cout << endl;

            cout << "friction:   " << fric << endl;
            cout << "young:    " << young << endl;
            cout << "poisson: " << poisson<< endl;
            cout << "Applied force " << force<< endl;
            cout << "Total  force step " << nstep<< endl;
            cout << "Total  volume step " << vnstep<< endl;
            cout << "Lateral parameter:f0 " << f0<< endl;
            cout << "spring constant " << kspr<< endl;
            cout << "spring radius gap " << sp_gap<< endl;

            cout << endl;
            cout << "visco_dissipation constant:   " << fric2 << endl;
            cout << "Angle dependancy: 0 (no depend), 1 (depend):  " << ang << endl;
            cout << "mesh refinement 1:yes " <<mesh_refine<< endl;

            cout << "Consistency check:  " << ConsCheck << endl;
            cout << "Friction control step:  " <<  vnstep+gap+control_fric<<" width: "<<width<< endl;

            cout << "Confinement potential:  " << Kconf<< endl;
            cout << "Total_moment_steps:  " << mom_steps<< endl;

            cout << "Surface tension:  " << gamma<< endl;
            cout << "Tension factor:  " << tens_factor<< "  tension plus :"  <<   gamma*2.0/(1.0+tens_factor)<<  " minus: " <<   gamma*2.0/(1.0+tens_factor) *tens_factor<<endl;

            cout << "Radius to thickness ratio:  " << R/height_in<< endl;

            cout << "tens_start:  " << tens_start<< endl;

            cout << "visco_elastic time:  " <<2*fric2*(1+ poisson)/young << endl;
            cout << "aap=1,bbn=1 normal mid:aap=2,bbn=0 up :aap=0,bbn=2 down::: " << aap << ": "<<bbn<< endl;
            cout << "width of transition zone:  " << width<< endl;
            cout << "bending constant:  " << kappa<< endl;

            cout << endl;
        }


    //bool  ConsCheck = false;

    // Time related parameters
    string sol_prefixMesh = "sol";
    double vol_step=1;
    RCP<DistributedMesh> posDisMesh,gloDisMesh,tensDisMesh;
    RCP<DOFsHandler> posDHand, gloDHand,viscoDHand,EDHand ;

    if (restart == 0)
    {

        // **************************************************************//
        // *****                   MESH CREATION                    *****//
        // **************************************************************//

        // Load mesh with MeshLoader
        RCP <MeshLoader> mesh = rcp(new MeshLoader);
        mesh->setMesh(ElemType::Triang, BasisFuncType::SubdivSurfs, 2);
        mesh->loadMesh(prefixMesh + ".vtk", MeshType::Sequential);



if(mesh_refine<1.0)
{
    // Distribute mesh for position
    posDisMesh = rcp(new DistributedMesh);
    posDisMesh->setMesh(mesh);
    posDisMesh->setBalanceMesh(true);
    posDisMesh->Update();
    posDisMesh->printFileLegacyVtk("pos_mesh_final_No_refine");

}
else
{

    // Distribute mesh for tension
    tensDisMesh = rcp(new DistributedMesh);
    tensDisMesh->setMesh(mesh);
    tensDisMesh->setBalanceMesh(true);
    tensDisMesh->Update();
    tensDisMesh->printFileLegacyVtk("pos_mesh_before_refine");


    // Distribute mesh for positions
    posDisMesh = rcp(new DistributedMesh);
    posDisMesh->setHRefinement(1);
    posDisMesh->setMeshRelation(MeshRelation::hRefin, tensDisMesh);
    posDisMesh->setBalanceMesh(true);
    posDisMesh->Update();
    posDisMesh->printFileLegacyVtk("pos_mesh_after_refined");

}

        // Distribute mesh global constraints
        gloDisMesh = rcp(new DistributedMesh);
        gloDisMesh->setMeshRelation(MeshRelation::GlobConstr, posDisMesh);
        gloDisMesh->setBalanceMesh(false);
        gloDisMesh->Update();

        gloDisMesh->printFileLegacyVtk("global_after");

        // **************************************************************//
        // *****               DOFHANDLER CREATION                  *****//
        // **************************************************************//
         posDHand = rcp(new DOFsHandler(posDisMesh));
        try
        {
            posDHand->setNameTag("posDHand");
            posDHand->setDOFs({"X", "Y", "Z"});
            posDHand->setNodeAuxF({"Ux", "Uy", "Uz","G11PA","G12PA","G22PA","G11NA","G12NA","G22NA","trGP","trGN","E1pE2","E1nE2","G1LA","G2LA","delWL_delC1","delWL_delC2","delWL_delG1","delWL_delG2","lam1","lam2","v11","v22","meanK","nx","ny","nz","height","jac"});
            posDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }



        //set initial condition;

        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);

            posDHand->nodeDOFs->setValue("X", i, IndexType::Local, x);
            posDHand->nodeDOFs->setValue("Y", i, IndexType::Local, y);
            posDHand->nodeDOFs->setValue("Z", i, IndexType::Local, z);

            posDHand->nodeAuxF->setValue(3, i, IndexType::Local, 1.0);
            posDHand->nodeAuxF->setValue(4, i, IndexType::Local, 0.0);
            posDHand->nodeAuxF->setValue(5, i, IndexType::Local, 1.0);
            posDHand->nodeAuxF->setValue(6, i, IndexType::Local, 1.0);
            posDHand->nodeAuxF->setValue(7, i, IndexType::Local, 0.0);
            posDHand->nodeAuxF->setValue(8, i, IndexType::Local, 1.0);
            posDHand->nodeAuxF->setValue(13, i, IndexType::Local, 1.0);//gl1
            posDHand->nodeAuxF->setValue(14, i, IndexType::Local, 1.0);//gl2
        }

        // InitializeNodalValues(posDHand, 3);
        posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);

        //Boundary conditions
        for(int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);

            int crease = posDisMesh->nodeCrease(i, hiperlife::IndexType::Local);

            if (spring<0.5)
            {

                if (choice ==10)
                {
                    if (sqrt(x*x+y*y) > R-sp_gap)
                    {
// boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }
                }
                else
                {
                    if (crease > 0)
                    {
// boundary condition

                        posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                        posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                    }
                }
            }

            if (out_choice > 0.1)
            {
                if (sqrt(x*x+y*y) > out_fact*R-sp_gap)
                {
// boundary condition

                    posDHand->setConstraint(0, i, hiperlife::IndexType::Local, 0.0);
                    posDHand->setConstraint(1, i, hiperlife::IndexType::Local, 0.0);
                    posDHand->setConstraint(2, i, hiperlife::IndexType::Local, 0.0);
                }
            }


        }



        posDHand->UpdateGhosts();
        posDHand->printFileLegacyVtk("checkinitialC", true);


// CREATE VISCO DOF HAND
        viscoDHand = rcp(new DOFsHandler(posDisMesh));
        try
        {
            viscoDHand->setNameTag("viscoDHand");
            viscoDHand->setDOFs({"G11P", "G12P", "G22P", "G11N", "G12N", "G22N", "G1L", "G2L"});
            viscoDHand->setNodeAuxF({"e1x", "e1y", "e1z", "e2x", "e2y", "e2z", "nx", "ny", "nz"});
            viscoDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": VISCO DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }

        LocMongeParam::computeLocal3DBasis(posDisMesh, viscoDHand->nodeAuxF);// assign auxiliary
        viscoDHand->UpdateGhosts();

        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            double x = posDisMesh->nodeCoord(i, 0, IndexType::Local);
            double y = posDisMesh->nodeCoord(i, 1, IndexType::Local);
            double z = posDisMesh->nodeCoord(i, 2, IndexType::Local);


            viscoDHand->nodeDOFs->setValue("G11P", i, IndexType::Local, 1.0);
            viscoDHand->nodeDOFs->setValue("G12P", i, IndexType::Local,0.0);
            viscoDHand->nodeDOFs->setValue("G22P", i, IndexType::Local, 1.0);
            viscoDHand->nodeDOFs->setValue("G11N", i, IndexType::Local,1.0);
            viscoDHand->nodeDOFs->setValue("G12N", i, IndexType::Local, 0.0);
            viscoDHand->nodeDOFs->setValue("G22N", i, IndexType::Local,1.0);

            viscoDHand->nodeDOFs->setValue("G1L", i, IndexType::Local, 1.0);
            viscoDHand->nodeDOFs->setValue("G2L", i, IndexType::Local,1.0);

        }
        // InitializeNodalValues(posDHand, 3);
        viscoDHand->nodeDOFs0->setValue(viscoDHand->nodeDOFs);

        viscoDHand->UpdateGhosts();
        viscoDHand->printFileLegacyVtk("checkinitialV", true);

        EDHand = rcp(new DOFsHandler(posDisMesh));
        try
        {
            EDHand->setNameTag("EDHand");
            EDHand->setDOFs({"EP", "EN","delWL_delC1","delWL_delC2","delWL_delG1","delWL_delG2","nx","ny","nz","height"});
            EDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }

        EDHand->UpdateGhosts();
        for (int i = 0; i < posDisMesh->loc_nPts(); i++)
        {
            EDHand->nodeDOFs->setValue("EP", i, IndexType::Local, 0.0);
            EDHand->nodeDOFs->setValue("EN", i, IndexType::Local,0.0);
            EDHand->nodeDOFs->setValue("delWL_delC1", i, IndexType::Local, 0.0);
            EDHand->nodeDOFs->setValue("delWL_delC2", i, IndexType::Local,0.0);
            EDHand->nodeDOFs->setValue("delWL_delG1", i, IndexType::Local, 0.0);
            EDHand->nodeDOFs->setValue("delWL_delG2", i, IndexType::Local,0.0);
        }

        EDHand->nodeDOFs0->setValue(EDHand->nodeDOFs);

        EDHand->UpdateGhosts();
        EDHand->printFileLegacyVtk("checkinitialE", true);



        gloDHand = rcp(new DOFsHandler(gloDisMesh));
        try
        {

            gloDHand->setNameTag("gloDHand");
            gloDHand->setDOFs({"P", "FX", "FY", "FZ", "MX", "MY", "MZ"});
            gloDHand->Update();
        }
        catch (runtime_error &err)
        {
            cout << myRank << ": DOFHandler could not be created " << err.what() << endl;
            MPI_Finalize();
            return 1;
        }

        gloDHand->setInitialCondition(0, 0.0);
        // gloDHand->setConstraint(0,0.0);
        gloDHand->setConstraint(1, 0.0);
        gloDHand->setConstraint(2, 0.0);
        gloDHand->setConstraint(3, 0.0);
        gloDHand->setConstraint(4, 0.0);
        gloDHand->setConstraint(5, 0.0);
        gloDHand->setConstraint(6, 0.0);
        gloDHand->UpdateGhosts();

        string fMesh = "sol_gloCons";
        gloDHand->printFile(fMesh, OutputMode::Text, true);


    }

    else
    {
        if (myRank == 0)
        cout<<"Starting from a restart" <<restart<<endl;
        //Cortex
        try
        {
            posDHand = rcp(new DOFsHandler("posDHand"));
            string fMesh = "sol_pos." + to_string(restart);
            posDHand->setFilePrefix(fMesh, OutputMode::Text);
            posDHand->setNumNodeAuxF(27);
            posDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":  Position DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }


//VISCO
        try
        {
            viscoDHand = rcp(new DOFsHandler("viscoDHand"));
            string fMesh = "sol_visco." + to_string(restart);
            viscoDHand->setFilePrefix(fMesh, OutputMode::Text);
            viscoDHand->setNumNodeAuxF(9);
            viscoDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":  visco DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }
//ENERGY HAND
        try
        {
            EDHand = rcp(new DOFsHandler("EDHand"));
            string fMesh = "sol_EN." + to_string(restart);
            EDHand->setFilePrefix(fMesh, OutputMode::Text);
            EDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":  ENERGY DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }



        //Global constraints
        try
        {
            gloDHand = rcp(new DOFsHandler("gloDHand"));
            string fMesh = "sol_gloCons." + to_string(restart);
            gloDHand->setFilePrefix(fMesh, OutputMode::Text);
            gloDHand->Update();
        }
        catch(runtime_error err)
        {
            cout << myRank << ":   gloCons DOFsHandler could not be created. " << err.what() << endl;

            MPI_Finalize();
            return 1;
        }


        //disMesh
        posDisMesh = posDHand->mesh;

    }


    //Resize auxiliary fields
    paramStr->a_aux.resize(posDisMesh->loc_nElem()*gPts*31);//just a definition
    paramStr->b_aux.resize(posDisMesh->loc_nElem()*gPts*8);



    // **************************************************************//
    // *****               HIPERPROBLEM CREATION                *****//
    // **************************************************************//
    SmartPtr<HiPerProblem> hiperProbl = Create<HiPerProblem>();
    try
    {
        // Set UserStructure
        hiperProbl->setParameterStructure(paramStr);
        hiperProbl->setConsistencyCheckDelta(1.E-6);
        hiperProbl->setConsistencyCheckTolerance(1.E-4);

        // Set DOFHandler
        hiperProbl->setDOFsHandlers({posDHand,gloDHand});

        // Set Integration
        hiperProbl->setIntegration("Integ", {"posDHand","gloDHand"});
        hiperProbl->setCubatureGauss("Integ", 3);
        //hiperProbl->setElementFillings("Integ", LS);
        if (ConsCheck == 1)
        {
            hiperProbl->setElementFillings("Integ", ConsistencyCheck<LS>);
        }
        else
        {
            hiperProbl->setElementFillings("Integ", LS);
        }
        hiperProbl->setGlobalIntegrals({"Energy","E1", "E2", "Dissipation", "Trd","volume","area","area_n","E3","Lat_area","Epress","E_tension1","E_tension2","E_tension3"});

        // Update
        hiperProbl->Update();
    }
    catch (runtime_error& err)
    {
        cout << myRank << ": HiPerProblem could not be created " << err.what() << endl;
        MPI_Finalize();
        return 1;
    }


    // create visco problem

    SmartPtr<HiPerProblem> viscoProbl= Create<HiPerProblem>();
    viscoProbl->setParameterStructure(paramStr);
    viscoProbl->setDOFsHandlers({viscoDHand});
    viscoProbl->setIntegration("Integ", {"viscoDHand"});
    viscoProbl->setCubatureGauss("Integ", 3);
    viscoProbl->setElementFillings("Integ", LS_ReacDif);
    viscoProbl->setGlobalIntegrals({"totBonds"});
    viscoProbl->Update();

    // create ENERGY problem

    SmartPtr<HiPerProblem> ENProbl= Create<HiPerProblem>();
    ENProbl->setParameterStructure(paramStr);
    ENProbl->setDOFsHandlers({EDHand});
    ENProbl->setIntegration("Integ", {"EDHand"});
    ENProbl->setCubatureGauss("Integ", 3);
    ENProbl->setElementFillings("Integ", LS_ED);
    ENProbl->Update();


    // **************************************************************//
    // *****                 SOLVERS' CREATION                  *****//
    // **************************************************************//


    RCP<MUMPSDirectLinearSolver> linSolver = rcp(new MUMPSDirectLinearSolver());
    linSolver->setHiPerProblem(hiperProbl);
    linSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
    linSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
    linSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
    linSolver->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
    linSolver->setDefaultParameters();
    linSolver->setWorkSpaceMemoryIncrease(1000);
    linSolver->Update();

    RCP<NewtonRaphsonNonlinearSolver> nonlinSolver = rcp(new NewtonRaphsonNonlinearSolver());
    nonlinSolver->setLinearSolver(linSolver);
    nonlinSolver->setMaxNumIterations(MAXITER_NR);
    nonlinSolver->setResTolerance(RESTOL_NR);
    nonlinSolver->setSolTolerance(SOLTOL_NR);
    nonlinSolver->setLineSearch(true);
    nonlinSolver->setConvRelTolerance(false);
    nonlinSolver->setPrintIntermInfo(true);
    nonlinSolver->setPrintSummary(false);
    nonlinSolver->setResMaximum(1E4);
    nonlinSolver->setSolMaximum(1E4);
    nonlinSolver->setExitRelMaximum(1E4);
    nonlinSolver->Update();

    RCP<MUMPSDirectLinearSolver>viscoDirSolver= rcp(new MUMPSDirectLinearSolver());
    viscoDirSolver->setHiPerProblem(viscoProbl);
    viscoDirSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
    viscoDirSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
    viscoDirSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
    viscoDirSolver->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
    viscoDirSolver->setDefaultParameters();
    viscoDirSolver->setWorkSpaceMemoryIncrease(1000);
    viscoDirSolver->Update();


    RCP<MUMPSDirectLinearSolver>ENDirSolver= rcp(new MUMPSDirectLinearSolver());
    ENDirSolver->setHiPerProblem(ENProbl);
    ENDirSolver->setMatrixType(MUMPSDirectLinearSolver::MatrixType::General);
    ENDirSolver->setAnalysisType(MUMPSDirectLinearSolver::AnalysisType::Parallel);
    ENDirSolver->setOrderingLibrary(MUMPSDirectLinearSolver::OrderingLibrary::Auto);
    ENDirSolver->setVerbosity(MUMPSDirectLinearSolver::Verbosity::None);
    ENDirSolver->setDefaultParameters();
    ENDirSolver->setWorkSpaceMemoryIncrease(1000);
    ENDirSolver->Update();


    hiperProbl->FillLinearSystem();


    viscoProbl->FillLinearSystem();
    ENProbl->FillLinearSystem();

    double a0{3.141};

    double  base_area= 3.141*R*R;

    v_target=v_target/hiperProbl->globalIntegral("area_n")*base_area;
    if (myRank == 0)
    {
        cout << " " << " Volume: " << hiperProbl->globalIntegral("volume") << endl;

        a0 = hiperProbl->globalIntegral("area_n");
        cout << "Target volume:  " << abs(v_target)*a0<< endl;

        cout << " " << " Area : phi: " << hiperProbl->globalIntegral("area_n")<< ": "<< hiperProbl->globalIntegral("area") << endl;
    }





    // globalintegral
    //Open file to write global integrals and write headers
    ofstream gIntegFile;
    gIntegFile.open ("globalIntegrals.dat");
    gIntegFile << "TS time deltat";
    //print
    for (auto g:  hiperProbl->globIntegralNames())
        gIntegFile << " " << g;
        gIntegFile << " pressure" ;
        gIntegFile << " Zmax" ;
        gIntegFile << " stress" ;
        gIntegFile << " strain" ;
    for (auto g:  viscoProbl->globalIntegralNames())
        gIntegFile << " " << g;

    gIntegFile << endl;

    double fstep=0;
   // double  factor=1;
    double astep=0;
    double mstep=0;
    double press0=0.0;
    double stress=0.0;
    double strain=0.0;
    double Z_max=1.0;
    // Time loop
    int timeStep=restart;

//print initial

    string solName = "sol_dis." + to_string(timeStep);
    posDHand->printFileLegacyVtk(solName,true);
    solName = "sol_visco." + to_string(timeStep);
    viscoDHand->printFileLegacyVtk(solName,true);

    double v0=v_target*a0;

    double deltatMax1=deltatMax;
     double gamma1=gamma;
    while ((tSimu < totalTime) and (timeStep < totalSteps))
    {
        // Print info
        if (myRank == 0)
        {
            cout<< "TS: " << timeStep + 1  << " Time " << tSimu << " of " << totalTime << " with deltat=" << deltat << endl;
            cout << "Starting Newton-Raphson iteration for membrane evolution" << endl;
        }


        // reduce volume in decreaments
        if(timeStep<vnstep)
        {
            factor=(timeStep+1.0)/vnstep*v_target;
            vstep=vstep+1;
            if (myRank == 0)
                cout<< "volume increased step: "<< vstep<<" fraction: "<<(timeStep+1.0)/vnstep<<endl;
            vol_inc=1;
        }

        if(timeStep>vnstep+gap-control_steps)
        {
            deltat *= stepFactor;
            deltatMax=control_dt;
            if (myRank == 0)
                cout << "time reduction step: " << timeStep-(vnstep+gap-control_steps)<<" dt: " << deltat<< " :max dt:  " << deltatMax <<endl;
        }

        if(timeStep>vnstep+gap+v_cont*v_rn_step)
        {
            deltat /= stepFactor;
            deltatMax=deltatMax1;

            if (myRank == 0)
                cout << "trelaxing pahse: " << timeStep<<" dt: " << deltat<< " :max dt:  " << deltatMax <<endl;
        }


        if(timeStep>vnstep+gap)
        {
            if (v_r_step<v_cont*v_rn_step)
            {
                factor=(1-(v_r_step+1.0)/v_rn_step)*v_target;
                factor=factor*exp(coeff*deltat*v_r_step);
                v_r_step=v_r_step+1;
                if (myRank == 0)
                    cout << "volume reduction step: " << v_r_step <<"by factor: "<<factor/v_target<< endl;

            }

            if(thin_shell>0.5) //to test thin-shell theory
            {
                paramStr->setRealParameter(MembParams::fric2,10000000.0*fric2);
            }

        }
        


        paramStr->setRealParameter(MembParams::factor,factor);
        paramStr->setRealParameter(MembParams::vol_inc,vol_inc);
        paramStr->setRealParameter(MembParams::timeStep,timeStep);




        if(timeStep>tens_start)
        {
            if(astep<aSteps)
            {

                paramStr->setRealParameter(MembParams::gamma,gamma1 * (astep+1.0)/aSteps);

                astep=astep+1;
            }
            if (myRank == 0)
                cout<< "Tension: "<< paramStr->getRealParameter(MembParams::gamma)<<endl;

        }

        if (deltat<0.000001)
        {
            if (myRank == 0)
                cout << "stopping time dt: " << endl;
            return 0;
        }


        // Prepare for solver
        posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
        gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
        hiperProbl->UpdateGhosts();
        //Prepare for solver
        viscoDHand->nodeDOFs0->setValue(viscoDHand->nodeDOFs);
        viscoDHand->UpdateGhosts();
        EDHand->nodeDOFs0->setValue(EDHand->nodeDOFs);
        EDHand->UpdateGhosts();

        // Solve the problem
        bool converged = nonlinSolver->solve();

        // Print info
        if (myRank == 0)
            cout << "Finished Newton-Raphson iteration" << endl;

        // Check convergence
        if (converged)
        {

            hiperProbl->FillLinearSystem();


            //Solve for visco problem
            if (myRank == 0)
                cout << "  Solving visco problem: "<<endl;
            bool viscoConverged = viscoDirSolver->solve();

            if(viscoConverged)
            {
                if (myRank == 0)
                cout << "visco Direct solver converged: " << endl;
            }
            else
            {
                if (myRank == 0)
                    cout << "visco Direct solver did not converge: " << endl;
            }
            // Update solutions
            viscoProbl->UpdateSolution();
            viscoProbl->FillLinearSystem();

            if (myRank == 0)
                cout << "  Solving Energy problem: "<<endl;
            // ENERGY PROBLEM
            bool ENConverged = ENDirSolver->solve();

            if(ENConverged)
            {
                if (myRank == 0)
                    cout << "ENERGY Direct solver converged: " << endl;
            }
            ENProbl->UpdateSolution();
            ENProbl->FillLinearSystem();


            tSimu += deltat;
            timeStep += 1;

            press0= gloDHand->nodeDOFs->getValue(0,0)/deltat;

            if(test_jac==1)
            {
                // find maximum z
                tensor<double, 1> Z_cord1(posDisMesh->loc_nPts());

                for (int i = 0; i < posDisMesh->loc_nPts(); i++)
                {
                    double z3 = posDHand->nodeDOFs->getValue(2, i, IndexType::Local);
                    Z_cord1(i) = z3;
                }


                //local maximum
                double Z_max1 = zz_max(Z_cord1, posDisMesh->loc_nPts());
                // 'global_max' will be the result on the root process.
                double global_max;
                //   MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

                // Perform the reduction operation MPI_MAX to find the global maximum.
                MPI_Reduce(&Z_max1, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                Z_max = global_max;
                if (myRank == 0)
                    cout << "The maximum Z level is: " << Z_max << endl;

                double rad_s = (R * R + Z_max * Z_max) / (2 * Z_max);
                stress = -0.5 * rad_s * press0;
                strain =  Z_max* Z_max/(R*R);
                if (myRank == 0)
                    cout << " " << " tension : " << stress << " areal strain : " << strain << endl;

            }


            for (int i = 0; i< posDHand->mesh->loc_nPts();i++)
            {

                double x0 = posDHand->nodeDOFs0->getValue(0,i,IndexType::Local);
                double y0 = posDHand->nodeDOFs0->getValue(1,i,IndexType::Local);
                double z0 = posDHand->nodeDOFs0->getValue(2,i,IndexType::Local);

                double x = posDHand->nodeDOFs->getValue(0,i,IndexType::Local);
                double y = posDHand->nodeDOFs->getValue(1,i,IndexType::Local);
                double z = posDHand->nodeDOFs->getValue(2,i,IndexType::Local);

                posDHand->nodeAuxF->setValue("Ux",i,IndexType::Local,(x-x0));
                posDHand->nodeAuxF->setValue("Uy",i,IndexType::Local,(y-y0));
                posDHand->nodeAuxF->setValue("Uz",i,IndexType::Local,(z-z0));
            //updated lagrangian
             /*   posDHand->mesh->_nodeData->setValue(0,i,IndexType::Local,x);
                posDHand->mesh->_nodeData->setValue(1,i,IndexType::Local,y);
                posDHand->mesh->_nodeData->setValue(2,i,IndexType::Local,z);*/
                double G11PR = viscoDHand->nodeDOFs->getValue(0,i,IndexType::Local);
                double G12PR = viscoDHand->nodeDOFs->getValue(1,i,IndexType::Local);
                double G22PR = viscoDHand->nodeDOFs->getValue(2,i,IndexType::Local);
                double G11NR = viscoDHand->nodeDOFs->getValue(3,i,IndexType::Local);
                double G12NR = viscoDHand->nodeDOFs->getValue(4,i,IndexType::Local);
                double G22NR = viscoDHand->nodeDOFs->getValue(5,i,IndexType::Local);

                double GL11 = viscoDHand->nodeDOFs->getValue(6,i,IndexType::Local);
                double GL22 = viscoDHand->nodeDOFs->getValue(7,i,IndexType::Local);

                posDHand->nodeAuxF->setValue("G11PA",i,IndexType::Local,G11PR);
                posDHand->nodeAuxF->setValue("G12PA",i,IndexType::Local,G12PR);
                posDHand->nodeAuxF->setValue("G22PA",i,IndexType::Local,G22PR);

                posDHand->nodeAuxF->setValue("G11NA",i,IndexType::Local,G11NR);
                posDHand->nodeAuxF->setValue("G12NA",i,IndexType::Local,G12NR);
                posDHand->nodeAuxF->setValue("G22NA",i,IndexType::Local,G22NR);

                posDHand->nodeAuxF->setValue("G1LA",i,IndexType::Local,GL11);
                posDHand->nodeAuxF->setValue("G2LA",i,IndexType::Local,GL22);

                posDHand->nodeAuxF->setValue("trGP",i,IndexType::Local,G11PR+G22PR);
                posDHand->nodeAuxF->setValue("trGN",i,IndexType::Local,G11NR+G22NR);



                double E1PE2 = EDHand->nodeDOFs->getValue(0,i,IndexType::Local);
                double E1NE2 = EDHand->nodeDOFs->getValue(1,i,IndexType::Local);
                double delWL_delC11 = EDHand->nodeDOFs->getValue(2,i,IndexType::Local);
                double delWL_delC22 = EDHand->nodeDOFs->getValue(3,i,IndexType::Local);
                double delWL_delG11 = EDHand->nodeDOFs->getValue(4,i,IndexType::Local);
                double delWL_delG22 = EDHand->nodeDOFs->getValue(5,i,IndexType::Local);
                posDHand->nodeAuxF->setValue("E1pE2",i,IndexType::Local,E1PE2);
                posDHand->nodeAuxF->setValue("E1nE2",i,IndexType::Local,E1NE2);
                posDHand->nodeAuxF->setValue("delWL_delC1",i,IndexType::Local,delWL_delC11);
                posDHand->nodeAuxF->setValue("delWL_delC2",i,IndexType::Local,delWL_delC22);
                posDHand->nodeAuxF->setValue("delWL_delG1",i,IndexType::Local,delWL_delG11);
                posDHand->nodeAuxF->setValue("delWL_delG2",i,IndexType::Local,delWL_delG22);
                //


                double nxx = EDHand->nodeDOFs->getValue(6,i,IndexType::Local);
                double nyy = EDHand->nodeDOFs->getValue(7,i,IndexType::Local);
                double nzz = EDHand->nodeDOFs->getValue(8,i,IndexType::Local);
                double height1 = EDHand->nodeDOFs->getValue(9,i,IndexType::Local);


                posDHand->nodeAuxF->setValue("nx",i,IndexType::Local,nxx);
                posDHand->nodeAuxF->setValue("ny",i,IndexType::Local,nyy);
                posDHand->nodeAuxF->setValue("nz",i,IndexType::Local,nzz);
                posDHand->nodeAuxF->setValue("height",i,IndexType::Local,height1);


            }

            posDHand->nodeDOFs0->setValue(posDHand->nodeDOFs);
            viscoDHand->nodeDOFs0->setValue(viscoDHand->nodeDOFs);
            EDHand->nodeDOFs0->setValue(EDHand->nodeDOFs);

            hiperProbl->UpdateGhosts();
            viscoProbl->UpdateGhosts();
            ENProbl->UpdateGhosts();

            if (timeStep%nPrint==0)
            {
                string solName;
                solName = sol_prefixMesh + to_string(timeStep);
                if (myRank == 0)
                    cout << "Printing file " <<  solName << endl;
                // ppHand->printFileLegacyVtk(solName,true);

            }

            if (nonlinSolver->numberOfIterations() < MAXITER_NR and deltat < deltatMax)
            {
                deltat /= stepFactor;
                if (myRank == 0)
                    cout << "dividing dt "  << endl;
            }

        }


        else
        {
            if (myRank == 0)
                std::cout << "ERROR::EXIT from Newton-Raphson but linear solver converged" << endl;

            // Modify time-step
            deltat *= stepFactor;

            // Restore initial values
            posDHand->nodeDOFs->setValue(posDHand->nodeDOFs0);
            gloDHand->nodeDOFs->setValue(gloDHand->nodeDOFs0);
            viscoDHand->nodeDOFs->setValue(viscoDHand->nodeDOFs0);
            EDHand->nodeDOFs->setValue(EDHand->nodeDOFs0);
        }



        // Compute factor to reduce the volume
        double v1 = hiperProbl->globalIntegral("volume");
        double a1 = hiperProbl->globalIntegral("area");

        if (myRank == 0)
        {
            cout << std::scientific << std::setprecision(6);
            cout << " " << " VOLUME: " << hiperProbl->globalIntegral("volume") << endl;
            cout << " " << " AREA jac: " << hiperProbl->globalIntegral("area") << endl;
            cout << " " << " Lateral area: " << hiperProbl->globalIntegral("Lat_area") << endl;
            cout << " " << " AREA_0 : jacR: " << hiperProbl->globalIntegral("area_n") << endl;
            cout << " " << " AREA_ jac_C*jacR: " << hiperProbl->globalIntegral("Trd") << endl;
            cout << " " << " stretch: " << hiperProbl->globalIntegral("area")/hiperProbl->globalIntegral("area_n") << endl;

            cout << " " << " volume fraction: " << v1/v0<< endl;
           // cout << " " << " volume needed to be reduced by: " << pow(Vf,timeStep)<< endl;
            cout << " " << " Dissipation: " << hiperProbl->globalIntegral("Dissipation") << endl;
            cout << " " << " Pressure work: " << -1*hiperProbl->globalIntegral("Epress") << endl;

            cout << " " << " Total strain energy: " << hiperProbl->globalIntegral("Energy") << endl;
            cout << " " << " Lateral energy: " << hiperProbl->globalIntegral("E3") << endl;
            double p0= gloDHand->nodeDOFs->getValue(0,0)/deltat;
            cout << " " << " pressure : " << -p0 << endl;
	        cout << " " << "Excess area: " << abs(hiperProbl->globalIntegral("area")-hiperProbl->globalIntegral("area_n"))/base_area<< endl;

            cout << " " << " Energy from tension top: " << hiperProbl->globalIntegral("E_tension1") << " bottom: "<<   hiperProbl->globalIntegral("E_tension2")   <<endl;
            cout << " " << " Energy from tension lateral: " << hiperProbl->globalIntegral("E_tension3")<<endl;
            if(bound>0.5)
            {
                cout << " " << " perimeter: " << hiperProbl->globalIntegral("BorderLength") << endl;
                cout << " " << " total moment REFERENCE: " << hiperProbl->globalIntegral("BorderLengthR") << endl;
                cout << " " << " Angle: " << hiperProbl->globalIntegral("angle1") << endl;
            }


        }





        // Print solution
        if (timeStep%nPrint==0)
        {
            string solName = "sol_dis." + to_string(timeStep);
            posDHand->printFileLegacyVtk(solName,true);
            //solName = "sol_gloCons." + to_string(timeStep);
            //gloDHand->printFileLegacyVtk(solName,true);
           // solName = "sol_visco." + to_string(timeStep);
           // viscoDHand->printFileLegacyVtk(solName,true);

        }




        if (timeStep % nPrint == 0)
        {
            //In the time loop print global integrals
            //Write global integrals
            if (myRank == 0)
            {
                gIntegFile << timeStep << " " << tSimu << " " << deltat;
                for (auto g:  hiperProbl->globIntegrals())
                    gIntegFile << " " << g;
                    gIntegFile << " " << press0;
                    gIntegFile << " " << Z_max;
                    gIntegFile << " " << stress;
                    gIntegFile << " " << strain;
                for (auto g:  viscoProbl->globIntegrals())
                    gIntegFile << " " << g;
                gIntegFile << endl;
            }
        }

        if(timeStep>4600)
        {
            if (timeStep%(hl_print) ==0)
            {

                string solName = "sol_pos." + to_string(timeStep);
                posDHand->printFile(solName, OutputMode::Text, true, tSimu);


                solName = "sol_gloCons." + to_string(timeStep);
                gloDHand->printFile(solName, OutputMode::Text, true, tSimu);

                solName = "sol_visco." + to_string(timeStep);
                viscoDHand->printFile(solName, OutputMode::Text, true, tSimu);

                solName = "sol_EN." + to_string(timeStep);
                EDHand->printFile(solName, OutputMode::Text, true, tSimu);
            }
        }





        // Update time-related quantities
        paramStr->setRealParameter(MembParams::deltat,deltat);

    }



    // **************************************************************//
    // *****                    FINALIZE                        *****//
    // **************************************************************//

    MPI_Finalize();
    return 0;
}
