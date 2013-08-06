within FCSysTest;
package Subregions
  extends Modelica.Icons.ExamplesPackage;

  model Subregion "Test a single subregion"
    extends FCSys.Subregions.Examples.Subregion(
      'inclC+'=true,
      'inclSO3-'=true,
      'incle-'=true,
      inclH2=true,
      inclH2O=true,
      inclN2=true,
      inclO2=true);
    // Note:  H+ is excluded to prevent reactions.
    // TODO:  Create a separate model to test reactions.
    // Currently, there are no assertions.  This model just checks that the
    // simulation runs.

  end Subregion;

  model Test2Subregions
    "Test two subregions with an initial pressure difference"
    extends FCSys.Subregions.Examples.Subregions(
      n_x=0,
      'inclC+'=true,
      'incle-'=true,
      'inclSO3-'=true,
      inclH2=true,
      inclH2O=true,
      inclN2=true,
      inclO2=true,
      environment(final analysis=true));
    // Note:  H+ is excluded to prevent reactions.

    output FCSys.Quantities.Amount S(stateSelect=StateSelect.never) =
      subregion1.graphite.'C+'.S + subregion2.graphite.'C+'.S + subregion1.ionomer.
      'SO3-'.S + subregion2.ionomer.'SO3-'.S + subregion1.graphite.'e-'.S +
      subregion2.graphite.'e-'.S + subregion1.gas.H2.S + subregion2.gas.H2.S +
      subregion1.gas.H2O.S + subregion2.gas.H2O.S + subregion1.gas.N2.S +
      subregion2.gas.N2.S + subregion1.gas.O2.S + subregion2.gas.O2.S
      "Total entropy";

  equation
    assert(der(S) >= 0, "Entropy may not decrease.");
    annotation (experiment(StopTime=30));
  end Test2Subregions;

end Subregions;
