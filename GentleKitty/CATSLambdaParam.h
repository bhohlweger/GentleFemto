#ifndef GENTLEKITTY_CATSLAMBDAPARAM_H_
#define GENTLEKITTY_CATSLAMBDAPARAM_H_

/// \file CATSLambdaParam.h
/// \brief Definition of lambda parameter calculator
/// \author A. Mathis
#include <vector>
#include <iostream>

/// \class Particle
/// This class incorporates the single particle quantities to compute the lambda
/// parameters
/// This include: Purity and primary/feed-down fraction(s)
class Particle {
 public:
  /// Constructor
  Particle()
      : fPurity(0),
        fPrimary(0),
        fFeedDown(0) {
  }
  ;

  /// Constructor
  /// \param pur Purity of the particle
  /// \param prim Primary fraction of the particle
  /// \param sec FeedDown fraction(s) of the particle
  Particle(const double pur, const double prim, std::vector<double> sec)
      : fPurity(pur),
        fPrimary(prim),
        fFeedDown(sec) {
  }
  ;

  /// Destructor
  ~Particle() = default;

  /// Set the purity
  /// \param pur Purity of the particle
  void SetPurity(double pur) {
    fPurity = pur;
  }

  /// Set the primary fraction
  /// \param prim Primary fraction of the particle
  void SetPrimaryFraction(double prim) {
    fPrimary = prim;
  }

  /// Set the feed-down fraction(s)
  /// \param sec FeedDown fraction(s) of the particle
  void SetFeedDownFraction(std::vector<double> sec) {
    fFeedDown = sec;
  }

  /// Get the purity
  /// \return Purity of the particle
  double GetPurity() const {
    return fPurity;
  }

  /// Get the primary fraction
  /// \return Primary fraction of the particle
  double GetPrimaryFraction() const {
    return fPrimary;
  }

  double GetTotalFeedDownFraction() const {
    double feed = 0.;
    for (unsigned int i = 0; i < fFeedDown.size(); ++i) {
      feed += fFeedDown[i];
    }
    return feed;
  }

  /// Get the number of feed-down contributions
  /// \return Number of feed-down contributions
  unsigned int GetNumberOfFeedDownContributions() const {
    return fFeedDown.size();
  }

  /// Get the i-th feed-down contribution
  /// \return i-th feed-down contribution
  double GetFeedDownFraction(int i) const {
    return fFeedDown[i];
  }

  /// Get all feed-down contributions
  /// \return All feed-down contributions
  std::vector<double> GetFeedDown() const {
    return fFeedDown;
  }

 private:
  double fPurity;                  ///< Purity of the particle
  double fPrimary;                 ///< Primary fraction of the particle
  std::vector<double> fFeedDown;  ///< Feed -down fraction(s) of the particle
};

/// \class CATSLambdaParam
/// This class computes the lambda parameters with the input of particle objects
class CATSLambdaParam {
 public:
  enum Type {
    Primary = 0,   ///< Genuine particle
    FeedDown = 1,  ///< Feed-down from resonances, etc.
    Fake = 2,      ///< Misidentified particles
  };

  /// Constructor
  CATSLambdaParam()
      : fParticles(),
        fIsSame(false) {
  }
  ;

  /// Constructor
  /// \param part1 Particle 1
  /// \param part2 Particle 2
  /// \param isSame Flag whether the two particles are of the same kind
  CATSLambdaParam(const Particle &part1, const Particle &part2,
                  const bool isSame = false);

  /// Destructor
  ~CATSLambdaParam() = default;

  /// Set the particles
  /// \param pur1 Purity of particle 1
  void SetParticle(const Particle &part) {
    fParticles.emplace_back(part);
  }

  /// Dump the output
  void PrintLambdaParams() const;

  /// Get the lambda parameter
  /// \param type Type of the contribution to be evaluated
  /// \return Contribution of the genuine pair (works for primary and fakes)
  double GetLambdaParam(const Type type) const;

  /// Get the lambda parameter
  /// \param type1 Type of the contribution of particle 1
  /// \param type2 Type of the contribution of particle 2
  /// \param sec1 FeedDown contribution of particle 1
  /// \param sec2 FeedDown contribution of particle 2
  /// \return Contribution of the genuine pair
  double GetLambdaParam(const Type type1, const Type type2, const int sec1 = 0,
                        const int sec2 = 0) const;

  /// Check whether the unput from the particles is same (sums up to one, only two particles, ...)
  /// \return Error code in case of problems, 0 in case everything is fine
  int SanityCheck() const;

 private:
  std::vector<Particle> fParticles;  ///< The involved particles
  bool fIsSame;  ///< Flag whether the two particles are of the same kind
};

inline int CATSLambdaParam::SanityCheck() const {
  // Sanity checks
  if (fParticles.size() != 2) {
    std::cerr
        << "ERROR LambdaParam: Less or more than two particles in the container \n";
    return -1;
  }
  for (unsigned int i = 0; i < fParticles.size(); ++i) {
    if (std::abs(
        fParticles[i].GetPrimaryFraction()
            + fParticles[i].GetTotalFeedDownFraction() - 1.f) > 0.01) {
      std::cerr << "ERROR LambdaParam: The contributions of particle " << i
                << " do not sum up to unity \n";
      std::cerr << "  Primary fraction " << fParticles[i].GetPrimaryFraction()
                << " Feed-down fraction "
                << fParticles[i].GetTotalFeedDownFraction() << "\n";
      return -2;
    }
  }
  return 0;
}

#endif  // GENTLEKITTY_CATSLAMBDAPARAM_H_
