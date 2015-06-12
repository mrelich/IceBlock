#ifndef Constants_hh
#define Constants_hh

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Define constants to be used in calculation of e-field in //
// standard SI units.                                       //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

const G4double m_c    = 2.99792458e8;   // spead of light in vacuum
const G4double m_n    = 1.78;           // index of refraction for ice                  
const G4double m_nAir = 1.000293;       // index of refraction for ice
const G4double m_pi   = 3.141592653589; // Approximation of pi
const G4double m_mu   = 4*m_pi*1e-7;    // permeability 
const G4double m_e    = 1.602e-19;      // electric charge

const G4double m_tolerance = 1e-15;

const G4double m_CGS_to_SI_q = 3.336e-10; // esu to Coulomb
const G4double m_CGS_to_SI_c = 1e-2;       // cm to m
#endif
