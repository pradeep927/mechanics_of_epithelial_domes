/*
    *******************************************************************************
    Copyright (c) 2017-2021 Universitat Politècnica de Catalunya
    Authors: Daniel Santos-Olivan, Alejandro Torres-Sanchez and Guillermo Vilanova
    Contributors:
    *******************************************************************************
    This file is part of hiperlife - High Performance Library for Finite Elements
    Project homepage: https://git.lacan.upc.edu/HPLFEgroup/hiperlifelib.git
    Distributed under the MIT software license, see the accompanying
    file LICENSE or http://www.opensource.org/licenses/mit-license.php.
    *******************************************************************************
*/


#ifndef AUXPATFOR_H
#define AUXPATFOR_H


#include <iostream>
#include <mpi.h>

#include "hl_DistributedMesh.h"
#include "hl_FillStructure.h"
//#include "hl_UserStructure.h"
#include "hl_DOFsHandler.h"
#include "hl_HiPerProblem.h"

#include "hl_Math.h"
#include "hl_Array.h"
#include "hl_ParamStructure.h"

using namespace std;
using namespace hiperlife;

using Teuchos::rcp;
using Teuchos::RCP;



struct MembParams
{
    enum RealParameters
    {
        R,
        thick,
        Ainc,
        aap,
        kappa ,
        bbn ,
        width ,
        bound ,
        bc_const,
        ang_old ,
        mesh_refine,
        mom_steps ,
        kspr,
        sp_gap ,
        a11,
        a12,
        a33,
        width_kap,
        kap1 ,
         mu,
         nu,
         fric,
        young,
        poisson,
        force,
        Kconf,
        gamma,
       f0,
        g_ratio,
        v_target,
        vstep,
        v_r_step,
        coeff,
        ang,
      tens_factor,
        out_choice,
        pn,
        deltat,
        totalTime,
        stepFactor,
        factor,
        fric_fact,
        height_in,
        fric2,
        v_rn_step,
        vol_inc,
        control_dt,
        v_cont,
        out_fact,
        thin_shell

        };



    enum IntParameters
    {
        spring,
        aSteps,
        control_fric,
        vnstep,
        tens_start,
        choice,
        control_steps,
        gap,
        test_jac,
        timeStep,
        fric_start
    };


    HL_PARAMETER_LIST DefaultValues
    {
                    {"R", 0.0001},
                    {"thick", 0.1},
                    {"Ainc", 0.0},
                    {"aap", 1.0},
                    {"kappa", 0.0},
                    {"bbn", 1.0},
                    {"width", 40.0},
                    {"bound", 0.0},
                    {"bc_const", 1.0},
                    {"ang_old",0.0},
                    {"mesh_refine",0.0},
                    {"mom_steps",100.0},
                    {"kspr",0.0},
                   {"sp_gap",0.008},
                    {"a11",0.51},
                    {"a12",0.49},
                    {"a33", 0.98},
                    {"width_kap", 40.0},
                    {"kap1", 1.0},
                    {"pn", 1.0},
                    {"mu", 10.0},
                     {"nu", 0.25},
                    {"fric", 0.01},
                    {"young", 25.0},
                    {"poisson", 0.25},
                    {"force", 0.0},
                    {"Kconf", 0.0},
                    {"gamma", 0.0},
                    {"f0", 0.9},
                    {"g_ratio", 0.8},
                    {"v_target", 0.667},
                    {"vstep", 0.0},
                    {"v_r_step", 0.0},
                    {"coeff", 0.0},
                    {"ang", 0.0},
                    {"tens_factor", 1.0},
                    {"out_choice", 0.0},
                    {"deltat", 0.0},
                   {"totalTime", 0.0},
                   {"stepFactor", 0.0},
                    {"factor", 1.0},
                   {"fric_fact", 1.0},
                  {"height_in", 0.05},
                  {"fric2", 0.015},
                  {"v_rn_step", 1000.0},
                 {"vol_inc", 1.0},
                    {"control_dt", 0.0001},
                    {"v_cont", 0.98},
                   {"out_fact", 1.2},
                  {"thin_shell", 0.0},
               {"spring", 0},
                   {"vnstep", 500},
                  {"tens_start", 2},
                 {"choice", 100},
                 {"control_steps", 400},
                 {"gap", 2000},
                {"test_jac", 0},
                 {"aSteps", 0},
               {"control_fric", 0},
               {"timeStep", 0},
                {"fric_start", 0}
        //

                };

};






void LS_ReacDif(hiperlife::FillStructure & fillStr);
void LS_ED(hiperlife::FillStructure & fillStr);

void LS(hiperlife::FillStructure& fillStr);

#endif


void RHS_Border(hiperlife::FillStructure & fillStr);




double inline ConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double pot{};
    if(x(2)<0)
    {
        pot = -Kconf/3.0 * x(2)*x(2)*x(2);
    }

    return pot;
}

double inline dConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double dpot{};

    if(x(2)<0)
    {
        dpot =-1*Kconf*x(2)*x(2);
    }



    return dpot;
}

double  inline ddConfinementPotential(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double ddpot{};
    if(x(2)<0)
    {
        ddpot = -2.0*Kconf*x(2);
    }

    return ddpot;
}

double inline zz_max(hiperlife::Tensor::tensor<double,1>& Z_cord1,double size_n)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    // find minimum

    double  Z_max = Z_cord1(0);
    // search num in inputArray from index 0 to elementCount-1
    for (int i = 0; i < size_n; i++)
    {
        if (Z_cord1(i) > Z_max)
        {
            Z_max = Z_cord1(i);
        }
    }

    return  Z_max;
}

double inline ConfinementPotential1(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    double pot{};
    pot = Kconf/2.0 * x(0)*x(0);
    return pot;
}

double inline dConfinementPotential1(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double dpot{};

    dpot =Kconf*x(0);

    return dpot;
}

double  inline ddConfinementPotential1(hiperlife::Tensor::tensor<double,1>& x, double Kconf)
{
    using namespace hiperlife;
    using namespace hiperlife::Tensor;
    double ddpot{};

    ddpot = Kconf;


    return ddpot;
}



