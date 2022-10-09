// Copyright (C) 2020 EDF
// All Rights Reserved
// This code is published under the GNU Lesser General Public License (GNU LGPL)
#define BOOST_TEST_MODULE testFastCDF
#define BOOST_TEST_DYN_LINK
#include <vector>
#include <memory>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/chrono.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random.hpp>
#include <Eigen/Dense>
#include "StOpt/cdf/fastCDF.h"



using namespace std;
using boost::timer::cpu_timer;
using boost::timer::auto_cpu_timer;
using namespace StOpt;
using namespace Eigen;

double accuracyEqual = 1e-10;

/// \brief Calculate CDF with naive method
/// \param p_x particules (sample)  size : (dimension, nbSim)
/// \param p_z  rectilinear points
/// \param p_y  estimate for each p_x point
ArrayXd CDFDirect(const ArrayXXd &p_x, const vector< shared_ptr<ArrayXd> >    &p_z, const ArrayXd &p_y)
{
    // store nbpt per dimension before
    ArrayXi nbPtZ(p_z.size());
    nbPtZ(0) = 1;
    for (size_t id = 0; id < p_z.size() - 1; ++id)
        nbPtZ(id + 1) = nbPtZ(id) * p_z[id]->size();
    // nb point
    int nbPtZT = nbPtZ(p_z.size() - 1) * p_z[p_z.size() - 1]->size();
    // store index
    ArrayXi index(p_z.size());
    // for return
    ArrayXd retCDF = ArrayXd::Zero(nbPtZT);
    for (int ip  = 0; ip < nbPtZT  ; ++ip)
    {
        // index
        int ipt = ip;
        for (int id = p_z.size() - 1; id  >= 0; --id)
        {
            index(id) = static_cast<int>(ipt / nbPtZ(id));
            ipt -= index(id) * nbPtZ(id);
        }
        for (int is = 0 ; is < p_x.cols(); ++is)
        {
            bool bAdd = true;
            for (int id = 0; id < p_x.rows(); ++id)
            {
                if (p_x(id, is) > (*p_z[id])(index(id)))
                {
                    bAdd = false;
                    break;
                }
            }
            if (bAdd)
            {
                retCDF(ip) += p_y(is);
            }
        }
    }
    return retCDF / p_y.size();
}





BOOST_AUTO_TEST_CASE(testFastCDF)
{
    boost::mt19937 generator;
    boost::normal_distribution<double> normalDistrib;
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<double> > normalRand(generator, normalDistrib);


    for (int ndim = 2; ndim <= 6; ++ndim)
    {

        int nbSim =  2000 ;

        for (int iter = 0; iter < 4; ++iter)
        {

            // number of evaluation points per dimension
            int nbEvalPerDim = static_cast<int>(pow(nbSim, 1. / ndim));

            ArrayXXd x(ndim, nbSim);
            vector< shared_ptr< ArrayXd > > z(ndim) ;  // evaluation points

            for (int id = 0; id < ndim; ++id)
            {
                for (int ip = 0; ip < nbSim; ++ip)
                {
                    x(id, ip) =  normalRand();
                }
                ArrayXd zDim = ArrayXd::LinSpaced(nbEvalPerDim, x.row(id).minCoeff() + 0.01, x.row(id).maxCoeff() - 0.01);
                z[id] = make_shared<ArrayXd>(zDim);
            }
            // size of window
            ArrayXd h = ArrayXd::Constant(ndim, 0.2);


            // calculate CDF
            ArrayXd y = ArrayXd::Constant(nbSim, 1.);

            cpu_timer  timeNaive;
            timeNaive.start();
            // direct
            ArrayXd CDFNaive  =  CDFDirect(x, z, y);
            timeNaive.stop();
            boost::chrono::duration<double> secondNaive = boost::chrono::nanoseconds(timeNaive.elapsed().user);


            cpu_timer  timeFast;
            timeFast.start();
            // fast
            ArrayXd CDFFastR = fastCDF(x, z, y);
            timeFast.stop();
            boost::chrono::duration<double> secondFast = boost::chrono::nanoseconds(timeFast.elapsed().user);

            cout <<  " ndim " << ndim << " fast " <<  secondFast.count() << " naive  " <<  secondNaive.count() << endl ;

            for (int is = 0; is < CDFNaive.size(); ++is)
                BOOST_CHECK_CLOSE(CDFNaive(is), CDFFastR(is), accuracyEqual);

            nbSim *= 2;

        }
    }
}
