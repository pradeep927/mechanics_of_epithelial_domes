

#include "Amesos.h"


/// hiperlife headers
#include "hl_FillStructure.h"
#include "hl_Geometry.h"
#include "hl_SurfLagrParam.h"
#include "hl_Tensor.h"
#include "hl_LinearSolver_Direct_Amesos2.h"
#include <hl_LinearSolver_Direct_MUMPS.h>
#include <fstream>

/// Header to auxiliary functions
#include "AuxViscousInextMembrane.h"

void LS(hiperlife::FillStructure& fillStr)
{

    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    //------------------------------------------------------------------
    // [2] Input variables

    //[1.1] Positions
    auto &subFill = (fillStr)["posDHand"];
    int nDim = subFill.nDim;
    int pDim = subFill.pDim;
    int eNN = subFill.eNN;
    int numDOFs = subFill.numDOFs;

    tensor<double, 2,false> nborDOFs0(subFill.nborDOFs0.data(), eNN, numDOFs);
    tensor<double, 2,false> nborDOFs(subFill.nborDOFs.data(), eNN, numDOFs);
    tensor<double, 2,false> nborCoords(subFill.nborCoords.data(), eNN, nDim);

    tensor<double, 1,false> bf(subFill.nborBFs(), eNN);
    tensor<double, 2,false> Dbf(subFill.nborBFsGrads(), eNN, pDim);
    tensor<double, 3,false> DDbf(subFill.nborBFsHess(), eNN, pDim, pDim);

    //[1.3] Global constraints
    auto &g_subFill = (fillStr)["gloDHand"];
    int g_numDOFs = g_subFill.numDOFs;

    tensor<double, 1,false> gDOFs(g_subFill.nborDOFs.data(), g_numDOFs);
    double pressure = gDOFs(0);
    tensor<double, 1> F = gDOFs(range(1, 3));
    tensor<double, 1> L = gDOFs(range(4, 6));

    //[1.4] Parameters
    double deltat = fillStr.getRealParameter(MembParams::deltat);
    double fric = fillStr.getRealParameter(MembParams::fric);
    double factor = fillStr.getRealParameter(MembParams::factor);
    double young = fillStr.getRealParameter(MembParams::young) * deltat;;
    double nu = fillStr.getRealParameter(MembParams::nu);
    double thick = fillStr.getRealParameter(MembParams::thick);
    double Kconf = fillStr.getRealParameter(MembParams::Kconf) * deltat;
    double gamma = fillStr.getRealParameter(MembParams::gamma) * deltat;
    int tens_start = fillStr.getIntParameter(MembParams::tens_start);
    int timeStep = fillStr.getIntParameter(MembParams::timeStep);//timestep

    double R    = fillStr.getRealParameter(MembParams::R);
    double width   = fillStr.getRealParameter(MembParams::width);
    double kappa = fillStr.getRealParameter(MembParams::kappa) * deltat;
    double aap    = fillStr.getRealParameter(MembParams::aap);
    double bbn    = fillStr.getRealParameter(MembParams::bbn);

    double a11= fillStr.getRealParameter(MembParams::a11);
    double a12= fillStr.getRealParameter(MembParams::a12);
    double a33= fillStr.getRealParameter(MembParams::a33);
    double kap1= fillStr.getRealParameter(MembParams::kap1);
    double width_kap   = fillStr.getRealParameter(MembParams::width_kap);


    double fric_fact = fillStr.getRealParameter(MembParams::fric_fact);
    double fric2 = fillStr.getRealParameter(MembParams::fric2);
    double height_in = fillStr.getRealParameter(MembParams::height_in);
    double vol_inc = fillStr.getRealParameter(MembParams::vol_inc);
    int fric_start = fillStr.getIntParameter(MembParams::fric_start);
    int control_fric = fillStr.getIntParameter(MembParams::control_fric);

    double Lagrangian{};
    double f0 = fillStr.getRealParameter(MembParams::f0);
    double g_ratio = fillStr.getRealParameter(MembParams::g_ratio);
    double kspr= fillStr.getRealParameter(MembParams::kspr) * deltat;
    int spring= fillStr.getIntParameter(MembParams::spring);

    double sp_gap= fillStr.getRealParameter(MembParams::sp_gap);

    double tens_factor   = fillStr.getRealParameter(MembParams::tens_factor);
    double gamma_plus=gamma*2.0/(1.0+tens_factor);
    double gamma_minus=gamma*2.0/(1.0+tens_factor)*tens_factor;
    //visco variables

    double *auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID * subFill.cubaInfo.iPts + subFill.kPt) * 8];
    double G11PG = auxiliary_b[0];
    double G12PG = auxiliary_b[1];
    double G22PG = auxiliary_b[2];
    double G11NG = auxiliary_b[3];
    double G12NG = auxiliary_b[4];
    double G22NG = auxiliary_b[5];
    double GL11 = auxiliary_b[6];
    double GL22 = auxiliary_b[7];

    //OUTPUTS
    tensor<double, 2,false> Bp(fillStr.Bk(0).data(), eNN, numDOFs);
    tensor<double, 1,false> Bg(fillStr.Bk(1).data(), g_numDOFs);

    tensor<double, 4,false> App(fillStr.Ak(0, 0).data(), eNN, numDOFs, eNN, numDOFs);
    tensor<double, 3,false> Apg(fillStr.Ak(0, 1).data(), eNN, numDOFs, g_numDOFs);
    tensor<double, 3,false> Agp(fillStr.Ak(1, 0).data(), g_numDOFs, eNN, numDOFs);


    //------------------------------------------------------------------
    // [2] Compute variables
    // [2.1] Compute geometry in the reference
    tensor<double, 1> xRv(nDim);
    tensor<double, 1> xRu(nDim);
    tensor<double, 1> xRuu(nDim);
    tensor<double, 1> xRuv(nDim);
    tensor<double, 1> xRvv(nDim);
    tensor<double, 1> normalR(nDim);
    tensor<double, 2> metricR(pDim, pDim);
    tensor<double, 2> curvatureR(pDim, pDim);
    tensor<double, 3> cristSymR(pDim, pDim, pDim);
    SurfLagrParam::ChristoffelSymbols(cristSymR, curvatureR, metricR, normalR, xRu, xRv, xRuu, xRvv, xRuv, eNN,
                                      nborCoords, Dbf, DDbf);
    tensor<double, 1> xR = bf * nborCoords;
    tensor<double, 2> imetricR = metricR.inv();
    double jacR = sqrt(metricR.det());



    //[2.1] Previous time-step
    tensor<double, 1> xu_n(nDim);
    tensor<double, 1> xv_n(nDim);
    tensor<double, 1> normal_n(nDim);
    tensor<double, 2> metric_n(pDim, pDim);
    SurfLagrParam::Normal(normal_n, metric_n, xu_n, xv_n, eNN, nborDOFs0, Dbf);

    tensor<double, 1> xuu_n(nDim);
    tensor<double, 1> xuv_n(nDim);
    tensor<double, 1> xvv_n(nDim);
    tensor<double, 2> curva_n(pDim, pDim);
    tensor<double, 3> christsym_n(pDim, pDim, pDim);
    SurfLagrParam::ChristoffelSymbols(christsym_n, curva_n, metric_n, normal_n, xu_n, xv_n, xuu_n, xvv_n, xuv_n, eNN,
                                      nborDOFs0, Dbf, DDbf);

    tensor<double, 1> x_n = bf * nborDOFs0;
    tensor<double, 2> imetric_n = metric_n.inv();


    double jac_n = sqrt(metric_n.det());
    double xnormal_n = x_n * normal_n;

    // UNIFORMLY DISTRIBUTED FORCE
    // UNIFORMLY DISTRIBUTED FORCE
    tensor<double, 1> Coords = bf * nborCoords;


    double rad=sqrt(xR(0)*xR(0)+xR(1)*xR(1));

    if(rad<R-sp_gap)
    {
        kspr=0.0;
    }

    if(kap1>0.5)
    {
        kappa=kappa*(a11-a12* tanh(width_kap * (a33*R - rad)));

    }
    else
    {
        kappa=kappa*1.0;
    }


    //[2.2] Current time-step
    tensor<double, 1> xu(nDim);
    tensor<double, 1> xv(nDim);
    tensor<double, 1> xuu(nDim);
    tensor<double, 1> xuv(nDim);
    tensor<double, 1> xvv(nDim);
    tensor<double, 1> normal(nDim);
    tensor<double, 2> metric(pDim, pDim);
    tensor<double, 2> curva(pDim, pDim);
    tensor<double, 3> christsym(pDim, pDim, pDim);
    SurfLagrParam::ChristoffelSymbols(christsym, curva, metric, normal, xu, xv, xuu, xvv, xuv, eNN, nborDOFs, Dbf,
                                      DDbf);


    tensor<double, 1> x = bf * nborDOFs;


    if (timeStep > fric_start+control_fric)
    {
        fric = fric * fric_fact + (fric * fric_fact - fric) * tanh(width * (thick - x(2)));
    }

    tensor<double, 2> imetric = metric.inv();
    double jac = sqrt(metric.det());
    double meancurva = product(imetric, curva, {{0, 0},{1, 1}});
    double xnormal = x * normal;
    double xRnormalR = xR * normalR;
    //3D unit tensor
    tensor<double, 2> Id(3, 3);
    Id(0, 0) = 1;
    Id(0, 1) = 0;
    Id(0, 2) = 0;
    Id(1, 0) = 0;
    Id(1, 1) = 1;
    Id(1, 2) = 0;
    Id(2, 0) = 0;
    Id(2, 1) = 0;
    Id(2, 2) = 1;

    // rate of deformation tensor
    tensor<double, 2> rodt = metric - metric_n;
    tensor<double, 2> rodt_CnCn = imetric_n * rodt * imetric_n;

    //[2.3] Other variables
    double tr_rodt = (jac - jac_n) / jac_n;


    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // [3.1] First derivatives

    tensor<double, 4> d_metric(eNN, numDOFs, pDim, pDim);
    SurfLagrParam::d_Metric(d_metric, Dbf, xu, xv);

    tensor<double, 4> d_metric_CC = product(product(d_metric, imetric, {{2, 0}}), imetric, {{2, 0}});
    tensor<double, 4> d_metric_CnCn = product(product(d_metric, imetric_n, {{2, 0}}), imetric_n, {{2, 0}});

    tensor<double, 4> d_imetric = -1.0 * d_metric_CC;

    tensor<double, 2> d_jac(eNN, numDOFs);
    SurfLagrParam::d_Jac(d_jac, Dbf, xu, xv, normal);


    tensor<double, 3> d_normal(eNN, numDOFs, nDim);
    SurfLagrParam::d_Normal(d_normal, Dbf, xu, xv, normal, jac, d_jac);

    tensor<double, 4> d_curva(eNN, numDOFs, pDim, pDim);
    SurfLagrParam::d_Curva(d_curva, DDbf, xuu, xuv, xvv, normal, d_normal);

    tensor<double, 2> d_meancurva = product(imetric, d_curva, {{0, 2},{1, 3}}) + product(curva, d_imetric, {{0, 2},{1, 3}});



    //------------------------------------------------------------------
    // [4] second derivatives

    tensor<double, 6> dd_metric(eNN, numDOFs, eNN, numDOFs, pDim, pDim);
    SurfLagrParam::dd_Metric(dd_metric, Dbf);

    tensor<double, 6> dd_metric_CC = product(product(dd_metric, imetric, {{4, 0}}), imetric, {{4, 0}});

    tensor<double, 6> dd_imetric = -1.0 * dd_metric_CC -product(product(d_metric, d_imetric, {{2, 2}}).transpose({0, 1, 3, 4, 5, 2}),imetric, {{5, 0}}) -product(product(d_metric, imetric, {{2, 0}}), d_imetric, {{2, 2}}).transpose({0, 1, 3, 4, 2, 5});

    tensor<double, 4> dd_jac(eNN, numDOFs, eNN, numDOFs);
    SurfLagrParam::dd_Jac(dd_jac, Dbf, xu, xv, normal, d_normal);

    tensor<double, 5> dd_normal(eNN, numDOFs, eNN, numDOFs, nDim);
    SurfLagrParam::dd_Normal(dd_normal, Dbf, jac, d_jac, dd_jac, normal, d_normal);

    tensor<double, 6> dd_curva(eNN, numDOFs, eNN, numDOFs, pDim, pDim);
    SurfLagrParam::dd_Curva(dd_curva, DDbf, xuu, xuv, xvv, d_normal, dd_normal);


    tensor<double, 4> dd_meancurva = product(imetric, dd_curva, {{0, 4},{1, 5}}) + product(d_imetric, d_curva, {{2, 2},{3, 3}}) +product(d_curva, d_imetric, {{2, 2},{3, 3}}) + product(curva, dd_imetric, {{0, 4},{1, 5}});



    // SHELL ANALYSIS

    //dphi0
    tensor<double, 2> Dphi_03(nDim, 2);
    Dphi_03(all, 0) = xRu(all);
    Dphi_03(all, 1) = xRv(all);
//dphi02.2
    tensor<double, 2> Dphi_0(2, 2);
    Dphi_0(all, 0) = xRu(range(0, 1));
    Dphi_0(all, 1) = xRv(range(0, 1));
// inverse of dphi0
    tensor<double, 2> iDphi_0(2, 2);
    iDphi_0 = Dphi_0.inv();

//dphi_t 3*2
    tensor<double, 2> Dx_uv(nDim, 2);
    Dx_uv(all, 0) = xu(all);
    Dx_uv(all, 1) = xv(all);
    tensor<double, 2> Dxn_uv(nDim, 2);
    Dxn_uv(all, 0) = xu_n(all);
    Dxn_uv(all, 1) = xv_n(all);

    tensor<double, 2> Dpsi(nDim, 2);
    Dpsi = product(Dx_uv, iDphi_0, {{1, 0}}); //F=dphi_t*dphi0_inv

    tensor<double, 2> Dpsi_n(nDim, 2);
    Dpsi_n = product(Dxn_uv, iDphi_0, {{1, 0}}); //F=dphi_t*dphi0_inv

    tensor<double, 2> DpsiT(2, 3);
    DpsiT = Dpsi.transpose({1, 0});

    tensor<double, 2> DpsiT_n(2, 3);
    DpsiT_n = Dpsi_n.transpose({1, 0});

    double jac0 = jacR;

    tensor<double, 2> metricC(pDim, pDim); //C_ab
    tensor<double, 2> metricC_n(pDim, pDim); //C_ab

    metricC = product(DpsiT, Dpsi, {{1, 0}});
    metricC_n = product(DpsiT_n, Dpsi_n, {{1, 0}});

    double jacC = sqrt(metricC.det()); //
    double jacC_n = sqrt(metricC_n.det()); //
    // cout<<"refrence:  " << jacC<< ": "  <<endl;
    double height0 = height_in / (jacC_n);//height_in/jac_n


    double pn = 1.0;

    // SWAP HERE AS G11PG, G11NG

    tensor<double, 2> GP(pDim, pDim);
    tensor<double, 2> GN(pDim, pDim);
    tensor<double, 2> GL(pDim, pDim);
    /*  G11PG=1.0;
      G12PG=0.0;
      G22PG=1.0;
      G11NG=1.0;
      G12NG=0.0;
      G22NG=1.0;*/

    GP(0, 0) = G11PG; // G ^AB // imetricR(0,0)
    GP(0, 1) = G12PG;
    GP(1, 0) = G12PG;
    GP(1, 1) = G22PG;
    GN(0, 0) = G11NG;
    GN(0, 1) = G12NG;
    GN(1, 0) = G12NG;
    GN(1, 1) = G22NG;

    GL(0, 0) = GL11;
    GL(0, 1) = 0.0;
    GL(1, 0) = 0.0;
    GL(1, 1) = GL22;

    tensor<double, 2> iGP = GP.inv(); //G_AB
    tensor<double, 2> iGN = GN.inv();
    tensor<double, 2> iGL = GL.inv();

    tensor<double, 2> metricP(pDim, pDim); //gp ^A_b
    tensor<double, 2> metricPR(pDim, pDim);//gp ^A_b
    tensor<double, 2> metricP0(pDim, pDim);//gp ^A_b

    tensor<double, 2> metricN(pDim, pDim); //gn ^a_b
    tensor<double, 2> metricNR(pDim, pDim);//gnr ^a_b
    tensor<double, 2> metricN0(pDim, pDim);//gn0 ^a_b

    metricP = product(iDphi_0, product((metric + aap* height0 * curva), iDphi_0, {{1, 0}}), {{0, 0}}); //Cp_ab
    metricPR = product(iDphi_0, product((metricR + aap* height0 * curvatureR), iDphi_0, {{1, 0}}), {{0, 0}}); //Cp_ab
    metricP0 = product(iDphi_0, product((metric_n +  aap*height0 * curva_n), iDphi_0, {{1, 0}}), {{0, 0}}); //Cp_ab


    metricN = product(iDphi_0, product((metric - bbn* height0 * curva), iDphi_0, {{1, 0}}), {{0, 0}}); //Cp_ab
    metricNR = product(iDphi_0, product((metricR -  bbn*height0 * curvatureR), iDphi_0, {{1, 0}}), {{0, 0}}); //Cp_ab
    metricN0 = product(iDphi_0, product((metric_n - bbn* height0 * curva_n), iDphi_0, {{1, 0}}), {{0, 0}}); //Cp_ab




    tensor<double,2> imetricP(pDim,pDim); //C+^ab
    tensor<double,2> imetricP0(pDim,pDim);
    tensor<double,2> imetricPR(pDim,pDim);

    tensor<double,2> imetricN(pDim,pDim); //C-^ab
    tensor<double,2> imetricN0(pDim,pDim);
    tensor<double,2> imetricNR(pDim,pDim);


    imetricP=metricP.inv();
    imetricP0=metricP0.inv();
    imetricN=metricN.inv();
    imetricN0=metricN0.inv();
    imetricPR=metricPR.inv();
    imetricNR=metricNR.inv();

    double jacP       = sqrt(metricP.det()); // J plus // sqrt(det gP)
    double jacPR       = sqrt(metricPR.det()); // J plus // sqrt(det gP)
    double jacP0       = sqrt(metricP0.det());


    double jacN       = sqrt(metricN.det()); // J plus // sqrt(det gP)
    double jacNR       = sqrt(metricNR.det()); // J plus // sqrt(det gP)
    double jacN0       = sqrt(metricN0.det());



    //REAL EVOLVING


    tensor<double,4> d_metricP(eNN,numDOFs,pDim,pDim);
    tensor<double,4> d_imetricP(eNN,numDOFs,pDim,pDim);
    tensor<double,6> dd_metricP(eNN,numDOFs,eNN,numDOFs,pDim,pDim);

    tensor<double,4> d_metricN(eNN,numDOFs,pDim,pDim);
    tensor<double,4> d_imetricN(eNN,numDOFs,pDim,pDim);
    tensor<double,6> dd_metricN(eNN,numDOFs,eNN,numDOFs,pDim,pDim);


    d_metricP= product(iDphi_0,product((d_metric+aap*height0*d_curva), iDphi_0,{{3,0}}) ,{{0,2}}).transpose({1,2,0,3})  ;    //M3_ab
    dd_metricP=product(iDphi_0,product((dd_metric+aap*height0*dd_curva), iDphi_0,{{5,0}}) ,{{0,4}}).transpose({1,2,3,4,0,5});// correct

    d_metricN=product(iDphi_0,product((d_metric-bbn*height0*d_curva), iDphi_0,{{3,0}}) ,{{0,2}}).transpose({1,2,0,3}) ; // correct
    dd_metricN=product(iDphi_0,product((dd_metric-bbn*height0*dd_curva), iDphi_0,{{5,0}}) ,{{0,4}}).transpose({1,2,3,4,0,5});// correct


    tensor<double,4> d_metricP_CC  = product(product(imetricP,d_metricP,{{1,2}}).transpose({1,2,0,3}),imetricP,{{3,0}});
    d_imetricP = -1.0 * d_metricP_CC;

    tensor<double,4> d_metricN_CC  = product(product(imetricN,d_metricN,{{1,2}}).transpose({1,2,0,3}),imetricN,{{3,0}});
    d_imetricN = -1.0 * d_metricN_CC;


//LATERAL

//LATERAL
    tensor<double,2> CL(pDim,pDim);
    CL(0,0)=jacC;
    CL(0,1)=0.0;
    CL(1,0)=0.0;
    CL(1,1)=1/(jacC*jacC);

    double jacCL=1/sqrt(jacC);
    double jacCL_n=1/sqrt(jacC_n);
    tensor<double,2> iCL=CL.inv();

// JACOBIAN AND DERIVATIVES
    //REAL EVOLVING

    tensor<double,2> d_jacP(eNN,numDOFs);
    d_jacP=0.5*jacP* product(imetricP,d_metricP,{{0,2},{1,3}});// Good
    tensor<double,2> d_jacN(eNN,numDOFs);
    d_jacN=0.5*jacN* product(imetricN,d_metricN,{{0,2},{1,3}});// Good

    tensor<double,2> d_jacC=d_jac/jacR;
    tensor<double,2> d_jacCL=-0.5/(pow(jacC,1.5))*d_jacC;

    //REAL EVOLVING
    tensor<double,4> dd_jacP(eNN,numDOFs,eNN,numDOFs);
    tensor<double,2> dd_jacP1(eNN,numDOFs);
    dd_jacP1=product(imetricP,d_metricP,{{0,2},{1,3}});// dgP/dx:gP^-1
    dd_jacP=0.5* jacP*product(imetricP,dd_metricP,{{0,4},{1,5}})+ 0.5*jacP*product(d_imetricP,d_metricP,{{2,2},{3,3}})+0.5*outer(d_jacP,dd_jacP1);// good correct

    tensor<double,4> dd_jacN(eNN,numDOFs,eNN,numDOFs);
    tensor<double,2> dd_jacN1(eNN,numDOFs);
    dd_jacN1=product(imetricN,d_metricN,{{0,2},{1,3}});// dgP/dx:gP^-1
    dd_jacN=0.5* jacN*product(imetricN,dd_metricN,{{0,4},{1,5}})+ 0.5*jacN*product(d_imetricN,d_metricN,{{2,2},{3,3}})+0.5*outer(d_jacN,dd_jacN1);// good correct

    tensor<double,4> dd_jacC=dd_jac/jacR;
    tensor<double,4> dd_jacCL=-0.5/(pow(jacC,1.5))*dd_jacC+0.75/(pow(jacC,2.5))*outer(d_jacC,d_jacC);


    // derivative of lateral strains
    tensor<double,4> d_CL(eNN,numDOFs,pDim,pDim);
    d_CL(all,all,0,0)=d_jacC;
    d_CL(all,all,1,0)=0.0*d_jacC;
    d_CL(all,all,0,1)=0.0*d_jacC;
    d_CL(all,all,1,1)=-2/(jacC*jacC*jacC)*d_jacC;

    tensor<double,6> dd_CL(eNN,numDOFs,eNN,numDOFs,pDim,pDim);
    dd_CL(all,all,all,all,0,0)=dd_jacC;
    dd_CL(all,all,all,all,0,1)=0.0*dd_jacC;
    dd_CL(all,all,all,all,1,0)=0.0*dd_jacC;
    dd_CL(all,all,all,all,1,1)=-2/(jacC*jacC*jacC)*dd_jacC+6/(jacC*jacC*jacC*jacC)*outer(d_jacC,d_jacC);


    //  viscoelastic metric

    double jacGP       = (GP.det()); // J of GPP (G11P*G22P-G12P*G12P)
    double jacGN       = (GN.det()); // J of GPP (G11P*G22P-G12P*G12P)
    double jacGL      =  (GL.det()); // J of GPP (G11P*G22P-G12P*G12P)
//first and third invariants
    double I1P=product(GP,metricP ,{{0,0},{1,1}}) ;
    double I1N=product(GN,metricN ,{{0,0},{1,1}}) ;
    double I1L=product(GL,CL ,{{0,0},{1,1}}) ;



//

    double I3P=jacP*jacP* jacGP;
    double I3N=jacN*jacN*jacGN;
    double I3L=jacCL*jacCL*jacGL;

    double I3P_SR=sqrt(jacGP)*jacP;
    double I3N_SR=sqrt(jacGN)*jacN;
    double I3L_SR=sqrt(jacGL)*jacCL;

    tensor<double,2> d_I1P(eNN,numDOFs);
    tensor<double,2> d_I1N(eNN,numDOFs);
    tensor<double,2> d_I3P(eNN,numDOFs);
    tensor<double,2> d_I3N(eNN,numDOFs);

    tensor<double,2> d_I3P_SR(eNN,numDOFs);
    tensor<double,2> d_I3N_SR(eNN,numDOFs);

    tensor<double,2> d_I1L(eNN,numDOFs);
    tensor<double,2> d_I3L(eNN,numDOFs);
    tensor<double,2> d_I3L_SR(eNN,numDOFs);

    d_I1P=product(GP,d_metricP ,{{0,2},{1,3}}) ;
    d_I1N=product(GN,d_metricN ,{{0,2},{1,3}}) ;
    d_I1L=product(GL,d_CL ,{{0,2},{1,3}}) ;

    d_I3P=2*jacP*d_jacP* jacGP;
    d_I3N=2*jacN*d_jacN* jacGN;
    d_I3L=2*jacCL*d_jacCL* jacGL;

    d_I3P_SR=sqrt(jacGP)*d_jacP;
    d_I3N_SR=sqrt(jacGN)*d_jacN;
    d_I3L_SR=sqrt(jacGL)*d_jacCL;


    tensor<double,4> dd_I1P(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I1N(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I1L(eNN,numDOFs,eNN,numDOFs);

    tensor<double,4> dd_I3P(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3N(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3L(eNN,numDOFs,eNN,numDOFs);

    tensor<double,4> dd_I3P_SR(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3N_SR(eNN,numDOFs,eNN,numDOFs);
    tensor<double,4> dd_I3L_SR(eNN,numDOFs,eNN,numDOFs);

    dd_I1P=product(GP,dd_metricP ,{{0,4},{1,5}}) ;
    dd_I1N=product(GN,dd_metricN ,{{0,4},{1,5}}) ;
    dd_I3P=(2*jacP*dd_jacP+2*outer(d_jacP,d_jacP))* jacGP;
    dd_I3N=(2*jacN*dd_jacN+2*outer(d_jacN,d_jacN))* jacGN;
    dd_I3P_SR=sqrt(jacGP)*dd_jacP;
    dd_I3N_SR=sqrt(jacGN)*dd_jacN;

    dd_I1L=product(GL,dd_CL ,{{0,4},{1,5}}) ;
    dd_I3L=(2*jacCL*dd_jacCL+2*outer(d_jacCL,d_jacCL))* jacGL;
    dd_I3L_SR=sqrt(jacGL)*dd_jacCL;


    double mu=0.5*young/(1+nu) ;
    double lambda=mu ;

    double jacPRR=jacP*jacR;
    double jacNRR=jacN*jacR;
    double jacLRR=jacCL*jacR;

    //VISCO-ELASTIC ENERGY
    double Evisco                     = (0.5*lambda*log(I3P_SR)*log(I3P_SR)-mu*log(I3P_SR)+0.5*mu*(I1P-2) )*jacPRR+ (0.5*lambda*log(I3N_SR)*log(I3N_SR)-mu*log(I3N_SR)+0.5*mu*(I1N-2) )*jacNRR;

    double Evisco1                     = (0.5*lambda*log(I3P_SR)*log(I3P_SR)-mu*log(I3P_SR)+0.5*mu*(I1P-2) )*jacPRR;
    double Evisco2                    = (0.5*lambda*log(I3N_SR)*log(I3N_SR)-mu*log(I3N_SR)+0.5*mu*(I1N-2) )*jacNRR;
    double Evisco3                    = (0.5*lambda*log(I3L_SR)*log(I3L_SR)-mu*log(I3L_SR)+0.5*mu*(I1L-2) )*f0*jacLRR;



    Bp(all,range(0,2)) += ((lambda*log(I3P_SR)-mu)/I3P_SR*d_I3P_SR+0.5*mu*d_I1P)*jacPRR;

    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3P_SR)-mu)/I3P_SR*dd_I3P_SR+0.5*mu*dd_I1P)*jacPRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3P_SR))+mu)/(I3P_SR*I3P_SR)*outer(d_I3P_SR,d_I3P_SR)*jacPRR;// term 1
    App(all,range(0,2),all,range(0,2))  +=outer(d_jacP,((lambda*log(I3P_SR)-mu)/I3P_SR*d_I3P_SR+0.5*mu*d_I1P) )*jacR;

    Bp(all,range(0,2)) +=((lambda*log(I3N_SR)-mu)/I3N_SR*d_I3N_SR+0.5*mu*d_I1N)*jacNRR;

    App(all,range(0,2),all,range(0,2))  += ((lambda*log(I3N_SR)-mu)/I3N_SR*dd_I3N_SR+0.5*mu*dd_I1N)*jacNRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3N_SR))+mu)/(I3N_SR*I3N_SR)*outer(d_I3N_SR,d_I3N_SR)*jacNRR;
    App(all,range(0,2),all,range(0,2))  +=   outer(d_jacN,((lambda*log(I3N_SR)-mu)/I3N_SR*d_I3N_SR+0.5*mu*d_I1N) )*jacR;//term 2
////////////////////////////////////////////////////

    Bp(all,range(0,2)) +=(0.5*lambda*log(I3P_SR)*log(I3P_SR)-mu*log(I3P_SR)+0.5*mu*(I1P-2) )*d_jacP*jacR;
    App(all,range(0,2),all,range(0,2))  +=(0.5*lambda*log(I3P_SR)*log(I3P_SR)-mu*log(I3P_SR)+0.5*mu*(I1P-2) )*dd_jacP*jacR;
    App(all,range(0,2),all,range(0,2))  +=outer(((lambda*log(I3P_SR)-mu)/I3P_SR*d_I3P_SR+0.5*mu*d_I1P),d_jacP)*jacR;

    Bp(all,range(0,2)) +=(0.5*lambda*log(I3N_SR)*log(I3N_SR)-mu*log(I3N_SR)+0.5*mu*(I1N-2))*d_jacN*jacR;

    App(all,range(0,2),all,range(0,2))  +=(0.5*lambda*log(I3N_SR)*log(I3N_SR)-mu*log(I3N_SR)+0.5*mu*(I1N-2))*dd_jacN*jacR;
    App(all,range(0,2),all,range(0,2))  +=outer(((lambda*log(I3N_SR)-mu)/I3N_SR*d_I3N_SR+0.5*mu*d_I1N),d_jacN )*jacR;


    //LATERAL
    Bp(all,range(0,2)) +=((lambda*log(I3L_SR)-mu)/I3L_SR*d_I3L_SR+0.5*mu*d_I1L)*f0*jacLRR;

    App(all,range(0,2),all,range(0,2))  +=((lambda*log(I3L_SR)-mu)/I3L_SR*dd_I3L_SR+0.5*mu*dd_I1L)*f0*jacLRR;
    App(all,range(0,2),all,range(0,2))  +=(lambda*(1-log(I3L_SR))+mu)/(I3L_SR*I3L_SR)*outer(d_I3L_SR,d_I3L_SR)*f0*jacLRR;
    App(all,range(0,2),all,range(0,2))  +=outer(d_jacCL,((lambda*log(I3L_SR)-mu)/I3L_SR*d_I3L_SR+0.5*mu*d_I1L))*f0*jacR;

    Bp(all,range(0,2)) +=(0.5*lambda*log(I3L_SR)*log(I3L_SR)-mu*log(I3L_SR)+0.5*mu*(I1L-2) )*f0*d_jacCL*jacR;

    App(all,range(0,2),all,range(0,2))  +=(0.5*lambda*log(I3L_SR)*log(I3L_SR)-mu*log(I3L_SR)+0.5*mu*(I1L-2))*f0*dd_jacCL*jacR;
    App(all,range(0,2),all,range(0,2))  +=outer(((lambda*log(I3L_SR)-mu)/I3L_SR*d_I3L_SR+0.5*mu*d_I1L),d_jacCL)*f0*jacR;


    //[6.2] Friction


    double Efric                        = 0.5*fric*jac_n*product(x-x_n,x-x_n,{{0,0}});
    Lagrangian+=Efric;


    Bp(all,range(0,2))                 += fric*jac_n*outer(bf,x-x_n);
    App(all,range(0,2),all,range(0,2)) += fric*jac_n*outer(outer(bf,Id),bf).transpose({0,1,3,2});


    //[6.3] Body force
    /*tensor<double,1> Id1(3);
    tensor<double,1> Id2(3);*/
    tensor<double,1> Id3(3);

   /* Id1(0)=1;
    Id1(1)=0;
    Id1(2)=0;

    Id2(0)=0;
    Id2(1)=1;
    Id2(2)=0;*/


    Id3(0)=0;
    Id3(1)=0;
    Id3(2)=1;


/*    tensor<double,2> d_x(eNN,numDOFs);
    tensor<double,2> d_y(eNN,numDOFs);*/
    tensor<double,2> d_z(eNN,numDOFs);


 /*   d_x=outer(bf,Id1);
    d_y=outer(bf,Id2);*/
    d_z=outer(bf,Id3);



    // bending

    
        Bp(all, range(0, 2)) += kappa*jacR * d_meancurva*meancurva;
        App(all, range(0, 2), all, range(0, 2)) += kappa*jacR *(dd_meancurva*meancurva+outer(d_meancurva,d_meancurva));


//spring potential
if(spring>0.5)
{
    Lagrangian += jacR *0.5*kspr*product(x-xR,x-xR,{{0,0}});


    Bp(all,0)                 += kspr*jacR*bf*(x(0)-xR(0));
    App(all,0,all,0) += kspr*jacR*outer(bf,bf);

    Bp(all,1)                 += kspr*jacR*bf*(x(1)-xR(1));
    App(all,1,all,1) += kspr*jacR*outer(bf,bf);

    Bp(all,2)                 += kspr*jacR*bf*(x(2)-xR(2));
    App(all,2,all,2) += kspr*jacR*outer(bf,bf);
}

    // [5.4]Adhesion Potential

    //Kconf=000000.0;
    Lagrangian += jac * ConfinementPotential(x, Kconf);

    Bp(all,range(0,2))                 += jac*d_z*dConfinementPotential(x, Kconf);
    Bp(all,range(0,2))        += d_jac*ConfinementPotential(x, Kconf);

    App(all,range(0,2),all,range(0,2)) +=outer(d_jac,d_z)*dConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=jac*outer(d_z,d_z)*ddConfinementPotential(x, Kconf);

    App(all,range(0,2),all,range(0,2)) +=dd_jac*ConfinementPotential(x, Kconf);
    App(all,range(0,2),all,range(0,2)) +=outer(d_z,d_jac)*dConfinementPotential(x, Kconf);



//surface tension


//surface tension
    double E_tens1=0.0;
    double E_tens2=0.0;
    double E_tens3=0.0;

    if(timeStep>tens_start)
    {
        Lagrangian += jacPR*jacR * gamma_plus+jacNR*jacR * gamma_minus+(gamma*f0)/sqrt(jacC)*jacR;

        E_tens1=(jacPR*jacR * gamma_plus)/deltat;
        E_tens2=(jacNR*jacR * gamma_minus)/deltat;
        E_tens3=(gamma*f0)/sqrt(jacC)*jacR/deltat;

        Bp(all, range(0, 2)) += d_jacP*jacR * gamma_plus+d_jacN*jacR * gamma_minus;
        App(all, range(0, 2), all, range(0, 2)) +=dd_jacP*jacR * gamma_plus+dd_jacN*jacR * gamma_minus;


        Bp(all, range(0, 2)) += (gamma*f0)*(-0.5)/pow(jacC_n,1.5)*d_jacC*jacR;
        App(all, range(0, 2), all, range(0, 2)) +=(gamma*f0)*(-0.5)/pow(jacC_n,1.5)*dd_jacC*jacR;


    }


    //[6.4] Volume constraint
    double Evol=pressure/3.0 *(jac * xnormal - factor*jac_n*xnormal_n);
    Lagrangian+= Evol;


    Bp += pressure / 3.0 * (d_jac * xnormal + jac * (d_normal * x) + jac * outer(bf, normal));


     Bg(0) += (jac * xnormal / 3) - factor * jacR * 1;
     App += pressure / 3.0 * (dd_jac * xnormal + outer(d_jac, d_normal * x) + outer(d_jac, outer(bf, normal)) +
                              outer(d_normal * x, d_jac) + outer(outer(bf, normal), d_jac) + jac * dd_normal * x +
                              jac * outer(d_normal, bf).transpose({0, 1, 3, 2}) +
                              jac * outer(bf, d_normal.transpose({2, 0, 1})));
     Apg(all, all, 0) += 1.0 / 3.0 * (d_jac * xnormal + jac * d_normal * x + jac * outer(bf, normal));




      Agp = Apg.transpose({2,0,1});
    //[7] Add contributions to global integrals
     fillStr.addToGlobalIntegral("E3",Evisco3/deltat);
     fillStr.addToGlobalIntegral("Dissipation",Efric/deltat);

     fillStr.addToGlobalIntegral("E_tension1",E_tens1);
    fillStr.addToGlobalIntegral("E_tension2",E_tens2);
    fillStr.addToGlobalIntegral("E_tension3",E_tens3);



    fillStr.addToGlobalIntegral("Energy",Evisco/deltat+Evisco3/deltat);
    fillStr.addToGlobalIntegral("E1",Evisco1/deltat);
    fillStr.addToGlobalIntegral("E2",Evisco2/deltat);
    fillStr.addToGlobalIntegral("Lat_area", f0*jacCL_n*jacR);
    fillStr.addToGlobalIntegral("Epress",Evol/deltat);
    // cout<< "king: " <<jacC<<endl;

    fillStr.addToGlobalIntegral("Trd",jacC*jacR);
    fillStr.addToGlobalIntegral("area_n", jacR);
    fillStr.addToGlobalIntegral("area", jac);
    fillStr.addToGlobalIntegral("volume", (1.0/3.0) * jac*xnormal);



    //This is to pass info to the other problem
    tensor<double,2> delI1P_delGP(pDim,pDim);//del I1 /del G ^bR
    tensor<double,2> delI1N_delGN(pDim,pDim);//
    tensor<double,2> delI1L_delGL(pDim,pDim);//
    tensor<double,2> delI1L_delCL(pDim,pDim);//

    delI1P_delGP=   (metricP)    ; //
    delI1N_delGN=   (metricN)    ; //_Br
    delI1L_delGL=CL;
    delI1L_delCL=GL;
    tensor<double,2> delI3P_delGP(pDim,pDim);//del I1 /del G ^bR
    tensor<double,2> delI3N_delGN(pDim,pDim);//
    tensor<double,2> delI3L_delGL(pDim,pDim);//
    tensor<double,2> delI3L_delCL(pDim,pDim);//

    delI3P_delGP=jacP*jacP*jacGP*iGP; //_Br
    delI3N_delGN=jacN*jacN*jacGN*iGN;
    delI3L_delGL=jacCL*jacCL*jacGL*iGL;
    delI3L_delCL=jacCL*jacCL*jacGL*iCL;

    tensor<double,2> delWP_delGP(pDim,pDim);//del j+ /del G ^bR
    tensor<double,2> delWN_delGN(pDim,pDim);//
    tensor<double,2> delWL_delGL(pDim,pDim);//
    tensor<double,2> delWL_delCL(pDim,pDim);//

    delWP_delGP=0.5*(lambda*log(I3P_SR)-mu)/I3P*delI3P_delGP+0.5*mu*delI1P_delGP;
    delWN_delGN=0.5*(lambda*log(I3N_SR)-mu)/I3N*delI3N_delGN+0.5*mu*delI1N_delGN;
    delWL_delGL=0.5*(lambda*log(I3L_SR)-mu)/I3L*delI3L_delGL+0.5*mu*delI1L_delGL;
    delWL_delCL=0.5*(lambda*log(I3L_SR)-mu)/I3L*delI3L_delCL+0.5*mu*delI1L_delCL;



    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*31];

    auxiliary_a[0]  = x(0);
    auxiliary_a[1]  = x(1);
    auxiliary_a[2]  = x(2);

    auxiliary_a[3]  =  metricP(0,0);
    auxiliary_a[4]  =  metricP(0,1);
    auxiliary_a[5]  =  metricP(1,0);
    auxiliary_a[6]  =  metricP(1,1);

    auxiliary_a[7]  =  metricN(0,0);
    auxiliary_a[8]  =  metricN(0,1);
    auxiliary_a[9]  =  metricN(1,0);
    auxiliary_a[10]  = metricN(1,1);


    auxiliary_a[11] = delWP_delGP(0,0);
    auxiliary_a[12] = delWP_delGP(0,1);
    auxiliary_a[13] = delWP_delGP(1,0);
    auxiliary_a[14] = delWP_delGP(1,1);

    auxiliary_a[15] = delWN_delGN(0,0);
    auxiliary_a[16] = delWN_delGN(0,1);
    auxiliary_a[17] = delWN_delGN(1,0);
    auxiliary_a[18] = delWN_delGN(1,1);



    auxiliary_a[19] = jacR;
    auxiliary_a[20] = jac;

    auxiliary_a[21] = Evisco1/deltat;
    auxiliary_a[22] =Evisco2/deltat;
    if (f0>0.001)
    {
        auxiliary_a[23] = 1 * delWL_delGL(0, 0);
        auxiliary_a[24] = 1 * delWL_delGL(1, 1);

        auxiliary_a[25] = 1 * delWL_delCL(0, 0) / deltat;
        auxiliary_a[26] = 1 * delWL_delCL(1, 1) / deltat;
    }
    else
    {
        auxiliary_a[23] = 0 * delWL_delGL(0, 0);
        auxiliary_a[24] = 0 * delWL_delGL(1, 1);

        auxiliary_a[25] = 0 * delWL_delCL(0, 0) / deltat;
        auxiliary_a[26] = 0 * delWL_delCL(1, 1) / deltat;
    }

    auxiliary_a[27] =normal(0);
    auxiliary_a[28] =normal(1);
    auxiliary_a[29] = normal(2);
    auxiliary_a[30] = height0;




}


void LS_ReacDif(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // --------------------------------------------
    // [0] Initialize
    // --------------------------------------------
    auto& subFill = (fillStr)["viscoDHand"];

    //Parameters
    int numDOFs = subFill.numDOFs;
    int nDim    = subFill.nDim;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;

    //Coordinates and degrees of freedom
    tensor<double,2,false> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    tensor<double,2,false> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    tensor<double,2,false> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    tensor<double,1,false>    bf(subFill.nborBFsDers(0),eNN);
    tensor<double,2,false>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
    tensor<double,3,false>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output
    tensor<double,2,false>  Bk(fillStr.Bk(0).data(),eNN,numDOFs);
    tensor<double,4,false>  Ak(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters
    // --------------------------------------------
    //[1.4] Parameters

    double fric_fact = fillStr.getRealParameter(MembParams::fric_fact);

    double deltat = fillStr.getRealParameter(MembParams::deltat);
    double fric    = fillStr.getRealParameter(MembParams::fric);
    double fric2    = fillStr.getRealParameter(MembParams::fric2);




    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*31];
    tensor<double,1,false>      x(&auxiliary_a[0],3);
    tensor<double,2,false> CP0(&auxiliary_a[3],2,2);
    tensor<double,2,false> CN0(&auxiliary_a[7],2,2);
    tensor<double,2,false> delWP_delGP0(&auxiliary_a[11],2,2);
    tensor<double,2,false> delWN_delGN0(&auxiliary_a[15],2,2);


    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];
    double E1_dt= auxiliary_a[21];
    double E2_dt= auxiliary_a[22];

    double delWL_delG1= auxiliary_a[23];
    double delWL_delG2= auxiliary_a[24];
    double delWL_delC1_dt= auxiliary_a[25];
    double delWL_delC2_dt= auxiliary_a[26];
    //SUPG
double G11P0{1.0}, G12P0{0.0},G22P0{1.0}, G11N0{1.0}, G12N0{0},G22N0{1.0} , G11P{1.0}, G12P{0.0},G22P{1.0}, G11N{1.0}, G12N{0},G22N{1.0}, G1L0{1.0},G2L0{1.0},G1L{1.0},G2L{1.0};



         G11P0  = nborDOFs0(all,0) * bf;
         G12P0  = nborDOFs0(all,1) * bf;
         G22P0  = nborDOFs0(all,2) * bf;
         G11N0  = nborDOFs0(all,3) * bf;
         G12N0  = nborDOFs0(all,4) * bf;
         G22N0  = nborDOFs0(all,5) * bf;

         G1L0  = nborDOFs0(all,6) * bf;
         G2L0  = nborDOFs0(all,7) * bf;

        G11P  = nborDOFs(all,0) * bf;
        G12P  = nborDOFs(all,1) * bf;
        G22P  = nborDOFs(all,2) * bf;
        G11N  = nborDOFs(all,3) * bf;
        G12N  = nborDOFs(all,4) * bf;
        G22N  = nborDOFs(all,5) * bf;

        G1L  = nborDOFs(all,6) * bf;
        G2L = nborDOFs(all,7) * bf;



    //Concentration values


//cout<<"step: "<<timestep<<"check: "<<G11P0<<endl;

    // FIXME: We need to remove this

    tensor<double,2> GP0(2,2);
    GP0(0,0)=G11P0;
    GP0(0,1)=G12P0;
    GP0(1,0)=G12P0;
    GP0(1,1)=G22P0;
    tensor<double,2> GN0(2,2);
    GN0(0,0)=G11N0;
    GN0(0,1)=G12N0;
    GN0(1,0)=G12N0;
    GN0(1,1)=G22N0;



    // Compute reactions
    double dt_eta=-1/fric2;

    double enP11=dt_eta*delWP_delGP0(0,0);
    double enP12=dt_eta*delWP_delGP0(0,1);
    double enP22=dt_eta*delWP_delGP0(1,1);

    double enN11=dt_eta*delWN_delGN0(0,0);
    double enN12=dt_eta*delWN_delGN0(0,1);
    double enN22=dt_eta*delWN_delGN0(1,1);

    double enL1=dt_eta*delWL_delG1;
    double enL2=dt_eta*delWL_delG2;

    //cout<<"check: "<<  " :"<<dt_eta*delWP_delGP0<<endl;


    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------

    //Rhs
    Bk(all,0) = jacR * bf * (G11P0 + enP11);
    Bk(all,1) = jacR * bf * (G12P0 + enP12);
    Bk(all,2) = jacR * bf * (G22P0 + enP22);
    Bk(all,3) = jacR * bf * (G11N0 + enN11);
    Bk(all,4) = jacR * bf * (G12N0 + enN12);
    Bk(all,5) = jacR * bf * (G22N0 + enN22);

    Bk(all,6) = jacR * bf * (G1L0 + enL1);
    Bk(all,7) = jacR * bf * (G2L0 + enL2);
    //Matrix
    Ak(all,0,all,0) = jacR * outer(bf, bf) ;//
    Ak(all,1,all,1) =jacR * outer(bf, bf);
    Ak(all,2,all,2) = jacR * outer(bf, bf);

    Ak(all,3,all,3) = jacR * outer(bf, bf);
    Ak(all,4,all,4) = jacR * outer(bf, bf);
    Ak(all,5,all,5) = jacR * outer(bf, bf);
    Ak(all,6,all,6) = jacR * outer(bf, bf);
    Ak(all,7,all,7) = jacR * outer(bf, bf);

    // --------------------------------------------
    // [4] Extra
    // --------------------------------------------
    //cout<<"print: "<<G1L<< ": " <<G2L<< endl;
    // This is to pass info to the other problem
    double* auxiliary_b = &fillStr.paramStr->b_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*8];
    auxiliary_b[0] = G11P;
    auxiliary_b[1] = G12P;
    auxiliary_b[2] = G22P;
    auxiliary_b[3] = G11N;
    auxiliary_b[4] = G12N;
    auxiliary_b[5] = G22N;
    auxiliary_b[6] = G1L;
    auxiliary_b[7] = G2L;
    // Global integrals
    fillStr.addToGlobalIntegral("totBonds",jacR);

}



void LS_ED(hiperlife::FillStructure& fillStr)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;

    // --------------------------------------------
    // [0] Initialize
    // --------------------------------------------
    auto& subFill = (fillStr)["EDHand"];

    //Parameters
    int numDOFs = subFill.numDOFs;
    int nDim    = subFill.nDim;
    int pDim    = subFill.pDim;
    int eNN     = subFill.eNN;
    double deltat = fillStr.getRealParameter(MembParams::deltat);

    //Coordinates and degrees of freedom
    tensor<double,2,false> nborCoords(subFill.nborCoords.data(),eNN,nDim);
    tensor<double,2,false> nborDOFs0(subFill.nborDOFs0.data(),eNN,numDOFs);
    tensor<double,2,false> nborDOFs(subFill.nborDOFs.data(),eNN,numDOFs);

    //Basis functions and derivatives
    tensor<double,1,false>    bf(subFill.nborBFsDers(0),eNN);
    tensor<double,2,false>   Dbf(subFill.nborBFsDers(1),eNN,pDim);
    tensor<double,3,false>  DDbf(subFill.nborBFsDers(2),eNN,pDim,pDim);

    //output
    tensor<double,2,false>  Bk1(fillStr.Bk(0).data(),eNN,numDOFs);
    tensor<double,4,false>  Ak1(fillStr.Ak(0,0).data(),eNN,numDOFs,eNN,numDOFs);


    // --------------------------------------------
    // [1] Model parameters

    // --------------------------------------------
    // [2] Auxiliary variables
    // --------------------------------------------
    //FIXME: this has to be loaded from some structure
    double* auxiliary_a = &fillStr.paramStr->a_aux[(subFill.loc_elemID*subFill.cubaInfo.iPts+subFill.kPt)*31];
    tensor<double,1,false>      x(&auxiliary_a[0],3);
    tensor<double,2,false> CP0(&auxiliary_a[3],2,2);
    tensor<double,2,false> CN0(&auxiliary_a[7],2,2);
    tensor<double,2,false> delWP_delGP0(&auxiliary_a[11],2,2);
    tensor<double,2,false> delWN_delGN0(&auxiliary_a[15],2,2);


    //Jacobian
    double jacR= auxiliary_a[19];
    double jac= auxiliary_a[20];
    double E1PP= auxiliary_a[21];
    double E1NN= auxiliary_a[22];

    double delWL_delG1= auxiliary_a[23]/deltat;
    double delWL_delG2= auxiliary_a[24]/deltat;
    double delWL_delC1= auxiliary_a[25];
    double delWL_delC2= auxiliary_a[26];

    double nxx= auxiliary_a[27];
    double nyy= auxiliary_a[28];
    double nzz= auxiliary_a[29];
    double height1= auxiliary_a[30];

    // --------------------------------------------
    // [3] Fill rhs and matrix
    // --------------------------------------------

//Rhs
    Bk1(all,0) = jacR * bf * (E1PP+E1NN);
    Bk1(all,1) = jacR * bf * (E1PP-E1NN);
    Bk1(all,2) = jacR * bf * (delWL_delC1);
    Bk1(all,3) = jacR * bf * (delWL_delC2);
    Bk1(all,4) = jacR * bf * (delWL_delG1);
    Bk1(all,5) = jacR * bf * (delWL_delG2);


    Bk1(all,6) = jacR * bf * (nxx);
    Bk1(all,7) = jacR * bf * (nyy);
    Bk1(all,8) = jacR * bf * (nzz);
    Bk1(all,9) = jacR * bf * (height1);

    //Matrix
    Ak1(all,0,all,0) = jacR * outer(bf, bf) ;//
    Ak1(all,1,all,1) =jacR * outer(bf, bf);
    Ak1(all,2,all,2) = jacR * outer(bf, bf) ;//
    Ak1(all,3,all,3) =jacR * outer(bf, bf);
    Ak1(all,4,all,4) = jacR * outer(bf, bf) ;//
    Ak1(all,5,all,5) =jacR * outer(bf, bf);
    Ak1(all,6,all,6) = jacR * outer(bf, bf) ;//
    Ak1(all,7,all,7) =jacR * outer(bf, bf);
    Ak1(all,8,all,8) = jacR * outer(bf, bf) ;//
    Ak1(all,9,all,9) =jacR * outer(bf, bf);


}


