#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "TestLambdaParam"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "CATSLambdaParam.h"

// Test1: All primary, full purity
BOOST_AUTO_TEST_CASE(Simple_test_1) {
  Particle p1(1.f, 1.f, { { 0 } });
  Particle p2(1.f, 1.f, { { 0 } });

  CATSLambdaParam params(p1, p2);
  const float primFrac = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(1.f, primFrac, 0.01);
  BOOST_CHECK_CLOSE(0.f, fakeFrac, 0.01);
}

// Test2: All primary, all fake
BOOST_AUTO_TEST_CASE(Simple_test_2) {
  Particle p1(0.f, 1.f, { { 0 } });
  Particle p2(0.f, 1.f, { { 0 } });

  CATSLambdaParam params(p1, p2);
  const float primFrac = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(0.f, primFrac, 0.01);
  BOOST_CHECK_CLOSE(1.f, fakeFrac, 0.01);
}

// Test3: Test secondaries
BOOST_AUTO_TEST_CASE(Simple_test_3) {
  Particle p1(1, 0.5, { { 0.3, 0.2 } });
  Particle p2(1., 1.f, { { } });

  CATSLambdaParam params(p1, p2);
  const float primFrac = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float secFrac1 = params.GetLambdaParam(CATSLambdaParam::FeedDown,
                                               CATSLambdaParam::Primary);
  const float secFrac2 = params.GetLambdaParam(CATSLambdaParam::FeedDown,
                                               CATSLambdaParam::Primary, 1);
  const float fakeFrac = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(0.5, primFrac, 0.01);
  BOOST_CHECK_CLOSE(0.3, secFrac1, 0.01);
  BOOST_CHECK_CLOSE(0.2, secFrac2, 0.01);
  BOOST_CHECK_CLOSE(0., fakeFrac, 0.01);
}

// Test4: Too few/many particles
BOOST_AUTO_TEST_CASE(Simple_test_4) {
  Particle p1(1.f, 0.5f, { { 0.25f, 0.25f } });

  CATSLambdaParam params;
  params.SetParticle(p1);
  const float primFrac1 = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac1 = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(-1.f, primFrac1, 0.01);
  BOOST_CHECK_CLOSE(-1.f, fakeFrac1, 0.01);

  Particle p2 = p1;
  Particle p3 = p1;

  params.SetParticle(p2);
  params.SetParticle(p3);
  const float primFrac2 = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac2 = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(-1.f, primFrac2, 0.01);
  BOOST_CHECK_CLOSE(-1.f, fakeFrac2, 0.01);
}

// Test5: Contributions don't sum up to unity
BOOST_AUTO_TEST_CASE(Simple_test_5) {
  Particle p1(1.f, 0.3f, { { 0.25f, 0.25f } });
  Particle p2(1.f, 0.3f, { { 0.25f, 0.25f } });

  CATSLambdaParam params(p1, p2);
  const float primFrac = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(-2.f, primFrac, 0.01);
  BOOST_CHECK_CLOSE(-2.f, fakeFrac, 0.01);
}

// Test6: Call with the wrong arguments
BOOST_AUTO_TEST_CASE(Simple_test_6) {
  Particle p1(0.f, 1.f, { { 0 } });
  Particle p2(0.f, 1.f, { { 0 } });

  CATSLambdaParam params(p1, p2);
  const float secFrac = params.GetLambdaParam(CATSLambdaParam::FeedDown);
  BOOST_CHECK_CLOSE(-4.f, secFrac, 0.01);
}

// p-p in p-Pb 5 TeV
BOOST_AUTO_TEST_CASE(pPb_p_p_5TeV) {
  double PurityProton = 0.984265;
  double PurityXi = 0.88;

  double protonPrim = 0.862814;
  double protonSec = 0.09603;

  double ProtonPrim = protonPrim;
  double arrayPercLamProton = protonSec / (1. - protonPrim);

  const unsigned NumChannels_p = 4;
  double* Purities_p = new double[NumChannels_p];
  double* Fraction_p = new double[NumChannels_p];

  Purities_p[0] = PurityProton;
  Purities_p[1] = PurityProton;
  Purities_p[2] = PurityProton;
  Purities_p[3] = 1. - PurityProton;

  Fraction_p[0] = ProtonPrim;
  Fraction_p[1] = (1. - ProtonPrim) * (arrayPercLamProton);
  Fraction_p[2] = (1. - ProtonPrim) * (1. - arrayPercLamProton);
  Fraction_p[3] = 1.;

  double lam_pp = Purities_p[0] * Fraction_p[0] * Purities_p[0] * Fraction_p[0];
  double lam_pp_pL = Purities_p[0] * Fraction_p[0] * Purities_p[1]
      * Fraction_p[1] * 2;
  double lam_pp_fake = Purities_p[3] * Purities_p[0]
      + Purities_p[0] * Purities_p[3] + Purities_p[3] * Purities_p[3];

  const Particle p1(
      PurityProton,
      protonPrim,
      { { (1. - ProtonPrim) * arrayPercLamProton, (1. - ProtonPrim)
          * (1. - arrayPercLamProton) } });
  const Particle p2 = p1;
  //
  CATSLambdaParam params(p1, p2, true);
  const float primFrac = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float secFrac = params.GetLambdaParam(CATSLambdaParam::Primary,
                                              CATSLambdaParam::FeedDown);
  const float fakeFrac = params.GetLambdaParam(CATSLambdaParam::Fake);

  BOOST_CHECK_CLOSE(lam_pp, primFrac, 0.01);
  BOOST_CHECK_CLOSE(lam_pp_pL, secFrac, 0.01);
  BOOST_CHECK_CLOSE(lam_pp_fake, fakeFrac, 0.01);
}

// p-Xi in p-Pb 5 TeV
BOOST_AUTO_TEST_CASE(pPb_p_Xi_5TeV) {
  double PurityProton = 0.984265;
  double PurityXi = 0.88;

  double protonPrim = 0.862814;
  double protonSec = 0.09603;

  double ProtonPrim = protonPrim;
  double arrayPercLamProton = protonSec / (1. - protonPrim);

  const unsigned NumChannels_p = 4;
  double* Purities_p = new double[NumChannels_p];
  double* Fraction_p = new double[NumChannels_p];

  Purities_p[0] = PurityProton;
  Purities_p[1] = PurityProton;
  Purities_p[2] = PurityProton;
  Purities_p[3] = 1. - PurityProton;

  Fraction_p[0] = ProtonPrim;
  Fraction_p[1] = (1. - ProtonPrim) * (arrayPercLamProton);
  Fraction_p[2] = (1. - ProtonPrim) * (1. - arrayPercLamProton);
  Fraction_p[3] = 1.;

  const double Xim1530_to_Xim = 0.32 * (1. / 3.);
  const double Xin1530_to_Xim = 0.32 * (2. / 3.);
  const double Omegam_to_Xim = 0.1;
  const double OmegamXim_BR = 0.086;

  const unsigned NumChannels_Xim = 5;
  double* Purities_Xim = new double[NumChannels_Xim];
  double* Fraction_Xim = new double[NumChannels_Xim];

  Purities_Xim[0] = PurityXi;
  Purities_Xim[1] = PurityXi;
  Purities_Xim[2] = PurityXi;
  Purities_Xim[3] = PurityXi;
  Purities_Xim[4] = 1. - PurityXi;

  Fraction_Xim[0] = 1. - Xim1530_to_Xim - Xin1530_to_Xim
      - Omegam_to_Xim * OmegamXim_BR;
  Fraction_Xim[1] = Xim1530_to_Xim;
  Fraction_Xim[2] = Xin1530_to_Xim;
  Fraction_Xim[3] = Omegam_to_Xim * OmegamXim_BR;
  Fraction_Xim[4] = 1.;

  const double lam_pXim = Purities_p[0] * Fraction_p[0] * Purities_Xim[0]
      * Fraction_Xim[0];

  const double lam_pXim_pXim1530 = Purities_p[0] * Fraction_p[0]
      * Purities_Xim[1] * Fraction_Xim[1];

  const double lam_pXim_fake = Purities_p[3] * Purities_Xim[0]
      + Purities_p[0] * Purities_Xim[4] + Purities_p[3] * Purities_Xim[4];

  const Particle proton(
      PurityProton,
      protonPrim,
      { { (1. - ProtonPrim) * arrayPercLamProton, (1. - ProtonPrim)
          * (1. - arrayPercLamProton) } });
  const Particle xi(
      PurityXi,
      1. - Xim1530_to_Xim - Xin1530_to_Xim - Omegam_to_Xim * OmegamXim_BR, { {
          Xim1530_to_Xim, Xin1530_to_Xim, Omegam_to_Xim * OmegamXim_BR } });

  CATSLambdaParam params(proton, xi);
  const float primFrac = params.GetLambdaParam(CATSLambdaParam::Primary);
  const float secFrac_pXim = params.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown);
  const float fakeFrac = params.GetLambdaParam(CATSLambdaParam::Fake);
  BOOST_CHECK_CLOSE(lam_pXim, primFrac, 0.01);
  BOOST_CHECK_CLOSE(lam_pXim_pXim1530, secFrac_pXim, 0.01);
  BOOST_CHECK_CLOSE(lam_pXim_fake, fakeFrac, 0.01);
}
