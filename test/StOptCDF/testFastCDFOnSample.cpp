// Copyright (C) 2020 EDF
// All Rights Reserved
// This code is published under the GNU Lesser General Public License (GNU LGPL)
#define BOOST_TEST_MODULE testFastCDFOnSample
#define BOOST_TEST_DYN_LINK
#include <vector>
#include <memory>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/chrono.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random.hpp>
#include <Eigen/Dense>
#include <fastCDFOnSample.h>



using namespace std;
using boost::timer::cpu_timer;
using boost::timer::auto_cpu_timer;
using namespace StOpt;
using namespace Eigen;

double accuracyEqual = 1e-10;


BOOST_AUTO_TEST_CASE(testFastCDFOnSample)
{
    boost::mt19937 generator;
    boost::normal_distribution<double> normalDistrib;
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<double> > normalRand(generator, normalDistrib);


    int ndimMin = 2;
    int ndimMax = 6;

    for (int ndim = ndimMin; ndim <= ndimMax; ++ndim)
    {
        int nbSim = 200;
        for (int id = 0; id < 4; ++id)
        {

            ArrayXXd ptr(ndim, nbSim);
            for (int is = 0; is < nbSim ; ++is)
            {
                for (int id = 0 ; id < ndim; ++id)
                    ptr(id, is) = normalRand();
            }

            // cout << " PTR " << ptr << endl ;
            cpu_timer  timeSlide;
            timeSlide.start();


            cpu_timer timerNaiveGrid;
            timerNaiveGrid.start();
            // CDF not renormalized
            ArrayXd  CDFVal(nbSim);
            for (int i = 0 ; i < nbSim ; ++i)
            {
                CDFVal(i) = 0;
                for (int j = 0; j < nbSim; ++j)
                {
                    if ((ptr.col(j) - ptr.col(i)).maxCoeff() <= 1e-12)
                    {
                        CDFVal(i) += 1;
                    }
                }
            }
            // normalize
            CDFVal /= nbSim;


            // calculate real CDF
            ArrayXd toAdd = ArrayXd::Constant(nbSim, 1.);

            ArrayXd cdfDivid = fastCDFOnSample(ptr, toAdd);


            timeSlide.stop();
            boost::chrono::duration<double> secondsSlide = boost::chrono::nanoseconds(timeSlide.elapsed().user);
            double timeSlideR  = secondsSlide.count() ;

            BOOST_CHECK_SMALL((cdfDivid - CDFVal).abs().maxCoeff(), accuracyEqual);
            cout << " ndim" << ndim << " Nb sim " << nbSim << " Slide " << timeSlideR << " ERROR " << (cdfDivid - CDFVal).abs().maxCoeff() <<   endl ;

            nbSim *= 2;
        }
    }
}
