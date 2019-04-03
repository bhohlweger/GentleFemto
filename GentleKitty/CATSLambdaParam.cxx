/// \file CATSLambdaParam.cxx
/// \brief Implementation of lambda parameter calculator
/// \author A. Mathis

#include "CATSLambdaParam.h"

CATSLambdaParam::CATSLambdaParam(const Particle &part1, const Particle &part2,
                                 const bool isSame)
    : fParticles(),
      fIsSame(isSame) {
  fParticles.emplace_back(part1);
  fParticles.emplace_back(part2);
}

void CATSLambdaParam::PrintLambdaParams() const {
  std::cout << "================= \n";
  std::cout << "LAMBDA PARAMETERS \n";
  std::cout << "================= \n\n";
  std::cout << "----------------- \n";
  printf("Primary fraction: %.2f\n\n", GetLambdaParam(Primary) * 100.f);

  const unsigned int sec1 = fParticles[0].GetNumberOfFeedDownContributions();
  const unsigned int sec2 = fParticles[1].GetNumberOfFeedDownContributions();

  std::cout << "----------------- \n";
  std::cout << "Number of feed-down fractions: \n";
  std::cout << "Particle 1: " << sec1 << "\n";
  std::cout << "Particle 2: " << sec2 << "\n\n";

  std::cout << "----------------- \n";
  std::cout << "FeedDown fractions: \n";
  // Primary 1 - FeedDown 2
  for (unsigned int i = 0; i < sec2; ++i) {
    printf("Prim1 - Sec2 (%i): %.2f\n", i,
           GetLambdaParam(Primary, FeedDown, 0, i) * 100.);
  }

  if (!fIsSame) {
    // FeedDown 1 - Primary 2
    for (unsigned int i = 0; i < sec1; ++i) {
      printf("Sec1 (%i) - Prim2: %.2f\n", i,
             GetLambdaParam(FeedDown, Primary, i, 0) * 100.);
    }
  }

  // Fake 1 - FeedDown 2
  for (unsigned int i = 0; i < sec2; ++i) {
    printf("Fake1 - Sec2 (%i): %.2f\n", i,
           GetLambdaParam(Fake, FeedDown, 0, i) * 100.);
  }

  if (!fIsSame) {
    // FeedDown 1 - Fake 2
    for (unsigned int i = 0; i < sec1; ++i) {
      printf("Sec1 (%i) - Fake2: %.2f\n", i,
             GetLambdaParam(FeedDown, Fake, i, 0) * 100.);
    }
  }

  // FeedDown 1 - FeedDown 2
  for (unsigned int i = 0; i < sec1; ++i) {
    int start = (fIsSame) ? i : 0;
    for (unsigned int j = start; j < sec2; ++j) {
      printf("Sec1 (%i) - Sec2 (%i): %.2f\n", i, j,
             GetLambdaParam(FeedDown, FeedDown, i, j) * 100.);
    }
  }

  std::cout << "\n----------------- \n";
  std::cout << "Fake fraction: " << GetLambdaParam(Fake) *100.f << "\n";
  std::cout << "Fake 1: " << GetLambdaParam(Fake, Primary, 0, 0) *100.f << "\n";
  std::cout << "Fake 2: " << GetLambdaParam(Primary, Fake, 0, 0) *100.f << "\n";
  std::cout << "Both  : " << GetLambdaParam(Fake, Fake, 0, 0) *100.f << "\n";
  std::cout << "----------------- \n";
}

double CATSLambdaParam::GetLambdaParam(const Type type) const {
  const int sane = SanityCheck();
  if (sane != 0) {
    return sane;
  }

  const double purity1 = fParticles[0].GetPurity();
  const double purity2 = fParticles[1].GetPurity();
  const double primary1 = fParticles[0].GetPrimaryFraction();
  const double primary2 = fParticles[1].GetPrimaryFraction();

  if (type == Primary) {
    return purity1 * primary1 * purity2 * primary2;
  } else if (type == FeedDown) {
    std::cerr
        << "ERROR LambdaParam: Specify the type of feed-down contribution \n";
    return -4.f;
  } else if (type == Fake) {
    return purity1 * (1.f - purity2) + (1.f - purity1) * purity2
        + (1.f - purity1) * (1.f - purity2);
  } else {
    std::cerr << "ERROR LambdaParam: Contribution type not defined \n";
    return -3.f;
  }
}

double CATSLambdaParam::GetLambdaParam(const Type type1, const Type type2,
                                       const int sec1, const int sec2) const {
  const int sane = SanityCheck();
  if (sane != 0) {
    return sane;
  }
  const double purity1 = fParticles[0].GetPurity();
  const double purity2 = fParticles[1].GetPurity();
  const double primary1 = fParticles[0].GetPrimaryFraction();
  const double primary2 = fParticles[1].GetPrimaryFraction();

  if (type1 == Primary && type2 == Primary) {
    // both are primaries
    return GetLambdaParam(Primary);
  } else if (type1 == Fake && type2 == Fake) {
    // all fake contributions
    return GetLambdaParam(Fake);
  } else if (type1 == Primary && type2 == Fake) {
    // First particle is primary, second is fake
    return purity1 * primary1 * (1.f - purity2);
  } else if (type1 == Fake && type2 == Primary) {
    // First particle is fake, second is primary
    return (1.f - purity1) * primary2 * purity2;
  } else if (type1 == Primary && type2 == FeedDown) {
    // First particle is primary, second is feed-down
    const double scaling = (fIsSame) ? 2. : 1.;
    return scaling * purity1 * primary1 * purity2
        * fParticles[1].GetFeedDownFraction(sec2);
  } else if (type1 == FeedDown && type2 == Primary) {
    // First particle is feed-down, second is primary
    const double scaling = (fIsSame) ? 2. : 1.;
    return scaling * purity1 * fParticles[0].GetFeedDownFraction(sec1) * purity2
        * primary2;
  } else if (type1 == Fake && type2 == FeedDown) {
    // First particle is fake, second is feed-down
    const double scaling = (fIsSame) ? 2. : 1.;
    return scaling * (1. - purity1) * primary1 * purity2
        * fParticles[1].GetFeedDownFraction(sec2);
  } else if (type1 == FeedDown && type2 == Fake) {
    // First particle is feed-down, second is fake
    const double scaling = (fIsSame) ? 2. : 1.;
    return scaling * purity1 * fParticles[0].GetFeedDownFraction(sec1)
        * (1. - purity2) * primary2;
  } else if (type1 == FeedDown && type2 == FeedDown) {
    // Both particles are feed-down
    const double scaling = (fIsSame && (sec1!=sec2)) ? 2. : 1.;
    return scaling * purity1 * fParticles[0].GetFeedDownFraction(sec1) * purity2
        * fParticles[1].GetFeedDownFraction(sec2);
  } else {
    std::cerr << "ERROR LambdaParam: Contribution type not defined \n";
    return -3.f;
  }
}
