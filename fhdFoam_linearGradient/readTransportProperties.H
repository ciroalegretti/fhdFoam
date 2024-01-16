singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta("beta", dimless/dimTemperature, laminarTransport);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

dimensionedScalar mu0("mu0", dimensionSet(1,1,-2,0,0,-2,0), laminarTransport);//     1.256e-6;		// vacuum mag. permeab.

dimensionedScalar rho("rho", dimensionSet(1,-3,0,0,0,0,0), laminarTransport);

// Pyromagnetic coefficient [1/K]
dimensionedScalar beta_m("beta_m", dimless/dimTemperature, laminarTransport);

dimensionedScalar chi0("chi0", dimless, laminarTransport);

dimensionedScalar Hmax("Hmax", dimensionSet(0,-1,0,0,0,1,0), laminarTransport);
