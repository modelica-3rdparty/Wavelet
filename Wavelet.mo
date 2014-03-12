within ;
package Wavelet "Modelica Wavelet Library"

  package Examples "Some examples to use this Wavelet Library"

  function displayWavelet
      "Display the wavelet available in this Wavelet Library"
    import Modelica.Constants.*;
    import Modelica.Utilities.Streams;
    import Modelica_LinearSystems2.Utilities.Plot;
    import Modelica_LinearSystems2.Math.Complex;

    input Records.wavletDefinition wp
        "Parameters for generating wavelet and scaling functions. See wavFunc() for detailed description of the parameters.";

    protected
    Records.wavFuncOut wOut; // output data of the function wavFunc()
    Real modulus[:];
    Real phase[:];
    Integer k;

  algorithm
    // get the data for display
      wOut:=Families.wavFunc(wp);

    //-----------------------------------------------------------------------------------
    // display the curves based on wavelet types
    //-----------------------------------------------------------------------------------
    // orthogonal wavelet, display one scaling and one wavelet functions
    if wOut.wavType == 1 or wOut.wavType==3 then
      Plot.diagram(
      Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=wOut.phi1, legend="Scaling function"),
                 Plot.Records.Curve(x=wOut.x, y=wOut.psi1, legend="Wavelet function")},
          heading=wOut.wavName));

    //-----------------------------------------------------------------------------------
    // biorthogonal wavelet, display two scaling and two wavelet functions
    elseif wOut.wavType == 2 then
      Plot.diagramVector({
        Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=wOut.phi1, legend="Scaling function"),
                 Plot.Records.Curve(x=wOut.x, y=wOut.psi1, legend="Wavelet function")},
          heading=wOut.wavName,
          xLabel="For decomposition (Nd="+String(wp.Nd)+")"),
        Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=wOut.phi2, legend="Scaling function"),
                 Plot.Records.Curve(x=wOut.x, y=wOut.psi2, legend="Wavelet function")},
          xLabel="For reconstruction (Nr=" +String(wp.Nr)+")")});

    //-----------------------------------------------------------------------------------
    // non-orthogonal and non-biorthogonal wavelets without scaling function, display only one wavelet function
    elseif wOut.wavType == 4 then
      Plot.diagram(
      Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=wOut.psi1, legend="Wavelet function")},
          heading=wOut.wavName));

    //-----------------------------------------------------------------------------------
    // complex wavelets without scaling function, display the real and imaginary parts of the wavelet function
    // with modular and phase curves
    elseif wOut.wavType == 5 then
      modulus := sqrt(wOut.psi1.*wOut.psi1 + wOut.psi2.*wOut.psi2);
      phase := fill(0, size(wOut.psi1,1));
      for k in 1:size(wOut.psi1,1) loop
        phase[k] := Modelica.Math.atan3(wOut.psi2[k], wOut.psi1[k], 0);
      end for;

      Plot.diagramVector({
        Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=wOut.psi1, legend="Real part"),
                 Plot.Records.Curve(x=wOut.x, y=wOut.psi2, legend="Imaginary part")},
          heading=wOut.wavName,
          xLabel="Wavelet function"),
        Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=modulus, legend="Modulus")},
          xLabel="Modulus"),
        Plot.Records.Diagram(
          curve={Plot.Records.Curve(x=wOut.x, y=phase, legend="Angle")},
          xLabel="Angle")});

    end if;

  annotation (Documentation(info="<html>
<p>
This example displays the wavelet and scaling functions of a wavelet family in one diagram. For a wavelet that has no scaling function,
only the wavelet function is displayed. For the complex wavelets, the real and imaginary parts are displayed with two curves.
The user is able to select a wavelet by its name and set its parameters through a graphic interface. 
</p><p>
Please refer to the description of
<a href=\"modelica://Wavelet.Families\">wavelet families</a> for the detailed information about the available wavelets
in this library.
</p><p>
Because of the length of discrete Meyer wavelet, the display of it takes a long time.
</p>
</html>"));
  end displayWavelet;

    function cwtChirpCurve "CWT of a chirp signal with outputs shown in curves"

        input Real[:] signal = {1.0000, 1.0000, 0.9999, 0.9996, 0.9987, 0.9968, 0.9933, 0.9877,
     0.9790, 0.9665, 0.9491, 0.9257, 0.8954, 0.8568, 0.8091, 0.7510,
     0.6818, 0.6008, 0.5075, 0.4020, 0.2847, 0.1566, 0.0194,-0.1245,
    -0.2720,-0.4192,-0.5614,-0.6932,-0.8089,-0.9022,-0.9672,-0.9981,
    -0.9901,-0.9397,-0.8452,-0.7074,-0.5297,-0.3186,-0.0837, 0.1623,
     0.4043, 0.6256, 0.8087, 0.9371, 0.9970, 0.9785, 0.8779, 0.6985,
     0.4518, 0.1572,-0.1588,-0.4646,-0.7266,-0.9127,-0.9970,-0.9634,
    -0.8096,-0.5491,-0.2111, 0.1616, 0.5174, 0.8026, 0.9702, 0.9879,
     0.8457, 0.5600, 0.1736,-0.2491,-0.6323,-0.9017,-1.0000,-0.9003,
    -0.6148,-0.1960, 0.2705, 0.6827, 0.9443, 0.9880, 0.7949, 0.4037,
    -0.0939,-0.5730,-0.9055,-0.9959,-0.8104,-0.3923, 0.1445, 0.6439,
     0.9522, 0.9671, 0.6743, 0.1591,-0.4140,-0.8520,-0.9988,-0.7935,
    -0.2997, 0.3091, 0.8072, 1.0000} "Vector, the signal to be analyzed";
        input Real[:] scales = {2,8,32} "Scales";

        input Records.wavletDefinition wIn
        "Input parameters for wavelet filters and functions";

    protected
        Real x0[:] = {i for i in 1:size(signal, 1)};
        Real[:,:] coefs "the transform coefficients";

        Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram[size(scales,1)+1]
        "the diagrams to be plotted";

    algorithm
        coefs := Transform.cwtn(signal, scales, wIn);

        // setup the diagrams, first the original signal
        diagram[1] := Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
            curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x0,
              y=signal,
              legend="Orignal Signal")},
            heading="CWT",
            xLabel="Point",
            yLabel="Amplitude");

        for i in 1:size(scales,1) loop
            // setup diagrams from scales = 1 to scales = max defined scales.
            diagram[i+1] := Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
              curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
                x=x0,
                y=coefs[i,:],
                legend="Scale = "+String(scales[i]))},
                xLabel="Point",
              yLabel="Amplitude");
        end for;

        // plot
        Modelica_LinearSystems2.Utilities.Plot.diagramVector(diagram=diagram);

    annotation (Documentation(info="<html>
<p>
This function carries out continuous wavelet transform of a chirp signal with graphic user interface. The data to be 
analyzed are given as default input values in the function. 
Bothe the original data and the transform result are shown with curves. 
</p><p>
Since simple default values are provided with this function, a direct execution of this 
function without any extra input parameters will deliver a complete MRA analysis for an example.
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed 
information about the available wavelets.
</p>
</html>"));
    end cwtChirpCurve;

    function cwtChirpImage
      "CWT of a chirp signal with outputs shown with curves and images"
      import Modelica_LinearSystems2.Utilities.Plot.Records.Diagram;
      import Modelica_LinearSystems2.Utilities.Plot.Records.Curve;

        input Real[:] signal = {1.0000, 1.0000, 0.9999, 0.9996, 0.9987, 0.9968, 0.9933, 0.9877,
     0.9790, 0.9665, 0.9491, 0.9257, 0.8954, 0.8568, 0.8091, 0.7510,
     0.6818, 0.6008, 0.5075, 0.4020, 0.2847, 0.1566, 0.0194,-0.1245,
    -0.2720,-0.4192,-0.5614,-0.6932,-0.8089,-0.9022,-0.9672,-0.9981,
    -0.9901,-0.9397,-0.8452,-0.7074,-0.5297,-0.3186,-0.0837, 0.1623,
     0.4043, 0.6256, 0.8087, 0.9371, 0.9970, 0.9785, 0.8779, 0.6985,
     0.4518, 0.1572,-0.1588,-0.4646,-0.7266,-0.9127,-0.9970,-0.9634,
    -0.8096,-0.5491,-0.2111, 0.1616, 0.5174, 0.8026, 0.9702, 0.9879,
     0.8457, 0.5600, 0.1736,-0.2491,-0.6323,-0.9017,-1.0000,-0.9003,
    -0.6148,-0.1960, 0.2705, 0.6827, 0.9443, 0.9880, 0.7949, 0.4037,
    -0.0939,-0.5730,-0.9055,-0.9959,-0.8104,-0.3923, 0.1445, 0.6439,
     0.9522, 0.9671, 0.6743, 0.1591,-0.4140,-0.8520,-0.9988,-0.7935,
    -0.2997, 0.3091, 0.8072, 1.0000} "Vector, the signal to be analyzed";

        input Real[:] scales = {i for i in 1:16} "Scales";
        input Records.wavletDefinition wIn
        "Input parameters for wavelet filters and functions";

    protected
        Real[size(scales,1),size(signal, 1)] coefs "the transform coefficients";
        Real[size(scales,1),size(signal, 1)] coefs1
        "the absolute value of transform coefficients";

        Real[size(scales,1),size(signal, 1)] x0;
        Real[size(scales,1),size(signal, 1)] y0;

    algorithm
        coefs := Transform.cwtn(signal, scales, wIn);

        for i in 1:size(scales,1) loop
            for j in 1:size(signal,1) loop
                x0[i,j] :=scales[i];
                y0[i,j] :=j;
            end for;
        end for;

        coefs1 := abs(coefs);

        //(x,y,z,nx,ny,nz,gx,gy,gz) := Plot3D.Utilities.SurfaceTest1(5);
        Plot3D.plotSurface({Plot3D.Records.ParametricSurfaceData(plotData=Plot3D.Records.PlotData(x=x0, y=y0, z=coefs1),
                                                                 styleData=Plot3D.Internal.StyleDataSurface(plotVectorField=false,
                                                                                                    style1=3))},
                           Plot3D.Records.CoordinateSystem(T=Plot3D.Internal.Transform(anglex=26,angley=-204,anglez=129.5,sy=2),
                                                           xAxisLabel="Scales",
                                                           yAxisLabel="Point",
                                                           zAxisLabel="Coefficients"),
                           true,
                           true,
                           1);

      // display the original signal
      Modelica_LinearSystems2.Utilities.Plot.diagram(Diagram(curve={Curve(x=(1:size(signal,1)), y=signal, legend="Signal")},xLabel="Point",yLabel="u",heading="Original signal"));

    annotation (Documentation(info="<html>
<p>
This function carries out continuous wavelet transform of a chirp signal with graphic user interface. The data to be 
analyzed are given as default input values in the function. 
The original data are shown with a curve. The transform result will be shown in an image with pseudo-color. 
For showing this image, the 'Plot3D' library has to be available.
</p><p>
Since simple default values are provided with this function, a direct execution of this 
function without any extra input parameters will deliver a complete MRA analysis for an example.
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed 
information about the available wavelets.
</p>

</html>"));
    end cwtChirpImage;

  function fileData_cwtn "CWT of the simulation result read from a MAT file"
      input String fileName = "testSignal2.mat"
        "The MAT file containing the signal to be analyzed"
                                                          annotation(Dialog(group="Input data"));
      input String signalName = "y" "The signal name to be analyzed" annotation(Dialog(group="Input data"));
      input Real t0=0.0 "Start time of the signal to be analyzed" annotation(Dialog(group="Input data"));
      input Real t1=1.0 "End time of the signal to be analyzed" annotation(Dialog(group="Input data"));
      input Real fs=500.0 "Sampling frequency of the signal for analysis. 
    This is used to generate an equidistant time grid with the expected time step. 
    If fs == 0, the simulation time grid will be used, even if it is none-equidistant.
    Since wavelet analysis has to be done in equidistant time grid, non-equidistant time grids
    usually cause incorrect results."  annotation(Dialog(group="Input data"));
      input Real[:] scales = {i for i in 1:4:64} "Scales" annotation(Dialog(group="Input data"));
      input Records.wavletDefinition wIn
        "Input parameters for wavelet filters and functions"
                                                           annotation(Dialog(group="Input data"));

    protected
      Integer n "the complete data length";
      Real t[1,:] "time axis of the simulation data";
      Real x[1,:] "the simulation result data";
      Integer t0id "index of t0 in the time grid";
      Integer t1id "index of t1 in the time grid";
      Real tSeg[:] "time grid of the interested data segment";
      Real xSeg[:] "the interested data segment";
      Integer nSeg "length of the data segment";
      Real[:,:] coefs "the transform coefficients";
      Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram
        "the diagrams to be plotted";
      Real[:,:] x0_3d;
      Real[:,:] y0_3d;

  algorithm
      // read simulation results
      n := readTrajectorySize(fileName);
      t := readTrajectory(fileName, {"Time"}, n);
      x := readTrajectory(fileName,{signalName},n);

      // get the original time grid between [t0, t1]
      t0id := General.findIndex(t[1,:], t0);
      t1id := General.findIndex(t[1,:], t1);

      // check the time range
      if t0id <= 0 then
        t0id := 1;
      end if;
      if t1id > n then
        t1id := n;
      end if;
      if t0id > n or t1id <= 0 then
        return;
      end if;

      // get the interested data
      if fs > 0 then
        tSeg := t[1,t0id]:1/fs:t[1,t1id];
        xSeg := General.interpL(t[1,:], x[1,:], tSeg);
      else
        tSeg := t[1,t0id : t1id];
        xSeg := x[1,t0id : t1id];
      end if;
      nSeg := size(xSeg,1);

      // cwt
      x0_3d := zeros(size(scales,1), nSeg);
      y0_3d := zeros(size(scales,1), nSeg);
      coefs := Transform.cwtn(xSeg, scales, wIn);

      // plot the 3-d display
      for i in 1:size(scales,1) loop
          for j in 1:nSeg loop
              x0_3d[i,j] :=scales[i];
              y0_3d[i,j] :=tSeg[j];
          end for;
      end for;

      Plot3D.plotSurface({Plot3D.Records.ParametricSurfaceData(plotData=Plot3D.Records.PlotData(x=x0_3d, y=y0_3d, z=abs(coefs)),
                                                                styleData=Plot3D.Internal.StyleDataSurface(plotVectorField=false,
                                                                                                  style1=3))},
                          Plot3D.Records.CoordinateSystem(T=Plot3D.Internal.Transform(anglex=26,angley=-204,anglez=129.5,sy=2),
                                                          xAxisLabel="Scales",
                                                          yAxisLabel="Time (s)",
                                                          zAxisLabel="Coefficients"),
                          true,
                          true,
                          1);

      // plot the original signal
      Modelica_LinearSystems2.Utilities.Plot.diagramVector(diagram={Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
            curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=tSeg,
              y=xSeg,
              legend="Orignal Signal")},
            heading=fileName,
            yLabel=signalName)});
  annotation (Documentation(info="<html>
<p>
This function carries out continuous wavelet transform with graphic user interface. The data to be 
analyzed are to be read from the *.MAT data file that is generated by Dymola. 
The original data are shown in a curve. The transform result is shown in an image with pseudo-color. 
For showing this image, the 'Plot3D' library has to be available.
</p><p>
Since simple default values are provided with this function, a direct execution of this 
function without any extra input parameters will deliver a complete MRA analysis for an example.
Due to large computational amount, this example will take some time (usually below 30 seconds)
to show the transform results.
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed 
information about the available wavelets.
</p>

</html>"));
  end fileData_cwtn;

  function fileDataMRA "MRA of the simulation result read from a MAT file"
      input String fileName = "testSignal2.mat"
        "The MAT file containing the signal to be analyzed"
                                                          annotation(Dialog(group="Input data"));
      input String signalName = "y" "The signal name to be analyzed" annotation(Dialog(group="Input data"));
      input Real t0=0.0 "Start time of the signal to be analyzed" annotation(Dialog(group="Input data"));
      input Real t1=1.0 "End time of the signal to be analyzed" annotation(Dialog(group="Input data"));
      input Real fs=200.0 "Sampling frequency of the signal for analysis. 
    This is used to generate an equidistant time grid with the expected time step. 
    If fs == 0, the simulation time grid will be used, even if it is none-equidistant.
    Since wavelet analysis has to be done in equidistant time grid, non-equidistant time grids
    usually cause incorrect results."  annotation(Dialog(group="Input data"));
    input Records.wavletDefinition wd
        "Wavelet definition. Note: Only (bi)orthogonal wavelets can be used for MRA. Valid parameters are wavID (<=7), Nd, Nr and range (dMeyer)"
      annotation(Dialog(group="Input data"));
    input Records.mraParameters mp "Parameters for MRA" annotation(Dialog(group="Input data"));

    protected
      Integer n "the complete data length";
      Real t[1,:] "time axis of the simulation data";
      Real x[1,:] "the simulation result data";
      Integer t0id "index of t0 in the time grid";
      Integer t1id "index of t1 in the time grid";
      Real tSeg[:] "time grid of the interested data segment";
      Real xSeg[:] "the interested data segment";

  algorithm
      // read simulation results
      n := readTrajectorySize(fileName);
      t := readTrajectory(fileName, {"Time"}, n);
      x := readTrajectory(fileName,{signalName},n);

      // get the original time grid between [t0, t1]
      t0id := General.findIndex(t[1,:], t0);
      t1id := General.findIndex(t[1,:], t1);

      // check the time range
      if t0id <= 0 then
        t0id := 1;
      end if;
      if t1id > n then
        t1id := n;
      end if;
      if t0id > n or t1id <= 0 then
        return;
      end if;

      // get the interested data
      if fs > 0 then
        tSeg := t[1,t0id]:1/fs:t[1,t1id];
        xSeg := General.interpL(t[1,:], x[1,:], tSeg);
      else
        tSeg := t[1,t0id : t1id];
        xSeg := x[1,t0id : t1id];
      end if;

    // call MRA GUI
    MRA.mraGUI(xSeg, wd, mp);

  annotation (Documentation(info="<html>
<p>
This function carries out the complete wavelet MRA with graphic user interface. The data to be 
analyzed are to be read from the *.MAT data file that is generated by Dymola. Please note that 
only the wavelets of types 1, 2 and 3 are possible for MRA. 
Otherwise running errors will occur. The results are displayed with multiple curves. 
It is possible to show either the wavelet coefficients at different 
levels or the signals that are reconstructed with the coefficients of different single levels. 
</p><p>
Since simple default values are provided with this function, a direct execution of this 
function without any extra input parameters will deliver a complete MRA analysis for an example.
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed 
information about the available wavelets.

</p>

</html>"));
  end fileDataMRA;

  function fileDataDenoising
      "Denoising of the simulation result read from a MAT file"
    import Modelica_LinearSystems2.Utilities.Plot.Records.Diagram;
    import Modelica_LinearSystems2.Utilities.Plot.Records.Curve;

        input String fileName = "testSignal2.mat"
        "The MAT file containing the signal to be analyzed"
                                                          annotation(Dialog(group="Input data"));
      input String signalName = "y" "The signal name to be analyzed" annotation(Dialog(group="Input data"));
      input Real t0=0.0 "Start time of the signal to be analyzed" annotation(Dialog(group="Input data"));
      input Real t1=1.0 "End time of the signal to be analyzed" annotation(Dialog(group="Input data"));
      input Real fs=200.0 "Sampling frequency of the signal for analysis. 
    This is used to generate an equidistant time grid with the expected time step. 
    If fs == 0, the simulation time grid will be used, even if it is none-equidistant.
    Since wavelet analysis has to be done in equidistant time grid, non-equidistant time grids
    usually cause incorrect results."  annotation(Dialog(group="Input data"));
    input Records.wavletDefinition wd
        "Wavelet definition. Note: Only (bi)orthogonal wavelets can be used for MRA. Valid parameters are wavID (<=7), Nd, Nr and range (dMeyer)"
      annotation(Dialog(group="Input data"));
    input Records.denoisParameters dp
        "Data and parameters for wavelet denoising"                               annotation(Dialog(group="Input data"));
    output Real y[:] "Denoised signal" annotation(Dialog(enable = false));
    output Real noise[:] "Noise removed from signal" annotation(Dialog(enable = false));
    output Real c[:]
        "Original wavelet coefficients of all levels, details see wavDec()"
                                                                          annotation(Dialog(enable = false));
    output Real cth[:]
        "Wavelet coefficients of all levels after applying thresholds"                annotation(Dialog(enable = false));
    output Integer len[:]
        "Lengths of coefficients for each levels, details see wavDec()"
                                                                       annotation(Dialog(enable = false));

    protected
      Integer n "the complete data length";
      Real t[1,:] "time axis of the simulation data";
      Real x[1,:] "the simulation result data";
      Integer t0id "index of t0 in the time grid";
      Integer t1id "index of t1 in the time grid";
      Real tSeg[:] "time grid of the interested data segment";
      Real xSeg[:] "the interested data segment";
  algorithm
      // read simulation results
      n := readTrajectorySize(fileName);
      t := readTrajectory(fileName, {"Time"}, n);
      x := readTrajectory(fileName,{signalName},n);

      // get the original time grid between [t0, t1]
      t0id := General.findIndex(t[1,:], t0);
      t1id := General.findIndex(t[1,:], t1);

      // check the time range
      if t0id <= 0 then
        t0id := 1;
      end if;
      if t1id > n then
        t1id := n;
      end if;
      if t0id > n or t1id <= 0 then
        return;
      end if;

      // get the interested data
      if fs > 0 then
        tSeg := t[1,t0id]:1/fs:t[1,t1id];
        xSeg := General.interpL(t[1,:], x[1,:], tSeg);
      else
        tSeg := t[1,t0id : t1id];
        xSeg := x[1,t0id : t1id];
      end if;

      // call denoising GUI
      Denoising.denoisGUI(xSeg, wd, dp);
  annotation (Documentation(info="<html>
<p>
This function carries out denoising operation to the one-dimensional data read from a specified *.MAT file, which is a simulation database 
generated by Dymola. 
The original data, the decomposed data in 
different levels, the denoised data in all levels and the denoised data are displayed with diagrams. 
In addition, the noise that is removed from the signal is displayed in a diagram, too.
</p><p>
The default input values provide an example to show the functionality of this function. 
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed 
information about the available wavelets.
</p>
</html>"));
  end fileDataDenoising;

    model testSignal1 "Generates a test signal for analysis"
        extends Modelica.Blocks.Interfaces.SignalSource;

    protected
        constant Real pi=Modelica.Constants.pi;
    equation
        y = 10*Modelica.Math.sin(2*pi*5*time) + Modelica.Math.sin(2*pi*20*time) + Modelica.Math.sin(2*pi*80*time);

    annotation (Documentation(info="<html>
<p>
This simple model generates a sum of three sinusoidal signals with different frequencies and magnitudes. 
</p><p>
After running this model in Dymola, all simulation data are saved in testSignal1.MAT file in the hard disk. 
</p>
</html>"));
    end testSignal1;

    model testSignal2 "Generates a test signal for analysis"
        extends Modelica.Blocks.Interfaces.SignalSource;

    protected
        constant Real pi=Modelica.Constants.pi;
    equation
        if time < 0.5 then
            y = Modelica.Math.sin(2*pi*10*time) + 0.2*Modelica.Math.sin(2*pi*50*time);
        else
            y = Modelica.Math.sin(2*pi*50*time) + 0.2*Modelica.Math.sin(2*pi*10*time);
        end if;

    annotation (Documentation(info="<html>
<p>
This simple model generates a sum of sinusoidal signals with different frequencies and magnitudes. 
</p><p>
After running this model in Dymola, all simulation data are saved in testSignal2.MAT file in the hard disk. 
</p>
</html>"));
    end testSignal2;

  annotation (Documentation(info="<html>
<p>
This section provides several examples. They serve as a starting point for the users to start their work with this library. 
Most examples can be accessed through a graphic interface for setting input parameters and displaying results.
</p><p>
To execute the examples that read simulation result data from MAT files, the two Modelica models (testSigna1 and testSignal2) 
have to be firstly run in Dymola, so that the MAT files are generated and stored in the hard disk.
</p>
</html>"));
  end Examples;

  package MRA "Application: One-dimensional multi-resolution analysis"

    function mraGUI
      "Graphic user interface for multi-resolution analysis (MRA)"
      import Modelica_LinearSystems2.Utilities.Plot.Records.Diagram;
      import Modelica_LinearSystems2.Utilities.Plot.Records.Curve;

      input Real u[:] = {0.9594,    3.6278,    7.6766,    8.3415,   11.8650,    9.9187,    7.3011,    2.1178,
       -0.6917,   -5.4113,   -8.2192,  -10.2508,  -10.9104,   -9.3952,   -5.9226,   -3.7118,
        2.6206,    5.5949,    9.1856,   10.1543,   10.1751,    6.4246,    3.8328,    0.5647,
       -2.8840,   -5.3850,   -9.2680,   -8.7256,   -8.9021,   -6.0738,   -3.8166,   -0.6568}
        "Signal to be analyzed"                                                                                      annotation(Dialog(group="Input data"));
      input Records.wavletDefinition wd
        "Wavelet definition. Note: Only (bi)orthogonal wavelets can be used for MRA. Valid parameters are wavID (<=7), Nd, Nr and range (dMeyer)"
        annotation(Dialog(group="Input data"));
      input Records.mraParameters mp "Parameters for MRA" annotation(Dialog(group="Input data"));

    protected
      Records.wavFuncOut wOut; // output data of the function wavFunc()
      Integer levels "decomposition levels, <= 12";
      Real rTune[:] "tuning ratio of all levels";
      Real c0[:]
        "original wavelet coefficients of all levels, details see wavDec()";
      Integer len[:]
        "lengths of coefficients for each levels, details see wavDec()";
      Real c1[:]
        "tuned wavelet coefficients of all levels, details see wavDec()";
      Integer nSig=size(u,1) "signal length";
      Diagram diagrams[12+2,2]
        "(N+2)x2 martix of diagrams. Note: Non-explicit definition of diagram dimensions does NOT work!";
      String str1 "temp string";
      String str2 "temp string";
      Integer k "temp integer";
      Integer j "temp integer";
      Real ra[:] "temp real array";

    algorithm
      // check the input parameters
      levels := mp.decLevel;
      if levels > 12 then
        levels :=12;
      end if;

      // get the wavelet filters
      wOut:=Families.wavFunc(wd);

      // --------------------------------------------------------------
      // decomposition
      (c0, len):=Transform.wavDec(u, levels, wOut.lod, wOut.hid); //wavelet multilevel decomposition

      // --------------------------------------------------------------
      // tune the coefficients
      rTune := {mp.rD1, mp.rD2, mp.rD3, mp.rD4, mp.rD5, mp.rD6, mp.rD7, mp.rD8, mp.rD9, mp.rD10, mp.rD11, mp.rD12};
      rTune := cat(1,rTune[1:levels], {mp.rA});
      c1 := tuneCoef(c0, len, rTune);

      // --------------------------------------------------------------
      // generate the diagram data
      // --------------------------------------------------------------
      // the original signal
      if mp.mraData==1 then
        str1 := "Original data and wavelet coefficients";
      else
        str1 := "Original and decomposited data";
      end if;
      diagrams[1,1] := Diagram(curve={Curve(x=(1:nSig), y=u, legend="Original data")},
                              xLabel="Point",yLabel="u",
                              heading=str1);

      // the reconstructed signal after coefficient tuning
      if mp.mraData==1 then
        str1 := "Reconstructed data and tuned coefficients";
      else
        str1 := "Reconstructed and tuned data";
      end if;
      diagrams[1,2] := Diagram(curve={Curve(x=(1:nSig), y=Transform.wavRec(c1, len, wOut.lor, wOut.hir), legend="Reconstructed data")},
                              xLabel="Point",yLabel="y",
                              heading=str1);

      // approximation data
      for k in 1:levels+1 loop  // loop for all levels, including approximation level
        for j in 1:2 loop // loop for original and tuned data
          // data type: coefficients or signals
          if mp.mraData == 1 then  // wavelet coefficients are to be displayed
            if j==1 then // original data
              ra := Transform.wavCoef1(c0, len, k);
              str2 := "Wavelet coefficients";
            else // tuned data
              ra := Transform.wavCoef1(c1, len, k);
              str2 := "Tuned coefficients";
            end if;

          else    // reconstructed signals are to be display
            if j==1 then // original data
              ra := Transform.wavRec1(c0, len, wOut.lor, wOut.hir, k);
              str2 := "Decomposed data";
            else  // tuned data
              ra := Transform.wavRec1(c1, len, wOut.lor, wOut.hir, k);
              str2 := "Reconstructed data with tuned coefficients";
            end if;

          end if;

          // for approximation level the label is different as that for detail levels
          if k==levels+1 then // approximation level
            str1 := "A" + String(levels);
          else  // detail levels
            str1 := "D" + String(k);
          end if;

          // generate the diagram data
          diagrams[levels-k+3,j] := Diagram(curve={Curve(x=(1:size(ra,1)), y=ra, legend=str2)},
                                     xLabel="Point",yLabel=str1);

        end for;
      end for;

      // --------------------------------------------------------------
      // display
      Modelica_LinearSystems2.Utilities.Plot.diagramMatrix(diagrams[1:levels+2,:]);

    annotation (Documentation(info="<html>
<p>
This function carries out the complete wavelet MRA with graphic user interface. Please note that 
only the wavelets of types 1, 2 and 3 are possible for MRA. 
Otherwise running errors will occur. The results are displayed with multiple curves. 
It is possible to show either the wavelet coefficients at different 
levels or the signals that are reconstructed with the coefficients of different single levels. 
</p><p>
Since simple default values are provided with this function, a direct execution of this 
function without any extra input parameters will deliver a complete MRA analysis for an example.
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed 
information about the available wavelets.

</p>

</html>"), Commands(executeCall=Wavelet.Examples.displayWavelet()));
    end mraGUI;

    function mra
      "Wavelet multi-resolution analysis (MRA)Graphic user interface for multi-resolution analysis (MRA)"
      import Modelica_LinearSystems2.Utilities.Plot.Records.Diagram;
      import Modelica_LinearSystems2.Utilities.Plot.Records.Curve;

      input Real u[:] = {0.9594,    3.6278,    7.6766,    8.3415,   11.8650,    9.9187,    7.3011,    2.1178,
       -0.6917,   -5.4113,   -8.2192,  -10.2508,  -10.9104,   -9.3952,   -5.9226,   -3.7118,
        2.6206,    5.5949,    9.1856,   10.1543,   10.1751,    6.4246,    3.8328,    0.5647,
       -2.8840,   -5.3850,   -9.2680,   -8.7256,   -8.9021,   -6.0738,   -3.8166,   -0.6568}
        "Signal to be analyzed"                                                                                      annotation(Dialog(group="Input data"));
      input Records.wavletDefinition wd
        "Wavelet definition. Note: Only (bi)orthogonal wavelets can be used for MRA. Valid parameters are wavID (<=7), Nd, Nr and range (dMeyer)"
        annotation(Dialog(group="Input data"));
      input Records.mraParameters mp "Parameters for MRA" annotation(Dialog(group="Input data"));
      output Real y[:,:] "Reconstructed signals with tuned coefficients in all levels. y[1,:] -- reconstructed signal, y[2,:] -- 1st detail level, 
    y[3,:] -- 2nd detail level, ..., y[N+1,:] -- N-th detail level, y[N+2,:] -- approximation level, where N is the number of 
    the decomposition levels.";

    protected
      Records.wavFuncOut wOut; // output data of the function wavFunc()
      Real y0[12+2,size(u,1)] "MRA result";
      Integer levels "decomposition levels, <= 12";
      Real rTune[:] "tuning ratio of all levels";
      Real c0[:]
        "original wavelet coefficients of all levels, details see wavDec()";
      Integer len[:]
        "lengths of coefficients for each levels, details see wavDec()";
      Real c1[:]
        "tuned wavelet coefficients of all levels, details see wavDec()";
      Integer k "temp integer";

    algorithm
      // check the input parameters
      levels := mp.decLevel;
      if levels > 12 then
        levels :=12;
      end if;

      // get the wavelet filters
      wOut:=Families.wavFunc(wd);

      // --------------------------------------------------------------
      // decomposition
      (c0, len):=Transform.wavDec(u, levels, wOut.lod, wOut.hid); //wavelet multilevel decomposition

      // --------------------------------------------------------------
      // tune the coefficients
      rTune := {mp.rD1, mp.rD2, mp.rD3, mp.rD4, mp.rD5, mp.rD6, mp.rD7, mp.rD8, mp.rD9, mp.rD10, mp.rD11, mp.rD12};
      rTune := cat(1,rTune[1:levels], {mp.rA});
      c1 := tuneCoef(c0, len, rTune);

      // --------------------------------------------------------------
      // reconstruct signals in all levels
      y0[1,:] := Transform.wavRec(c1, len, wOut.lor, wOut.hir);
      for k in 1:levels+1 loop  // loop for all levels
          y0[k+1,:] := Transform.wavRec1(c1, len, wOut.lor, wOut.hir, k);
      end for;

      // output
      y := y0[1:levels+2, :];

    annotation (Documentation(info="<html>
<p>
This function carries out the complete wavelet MRA and outputs decomposed signals in all levels after coefficient tuning. 
Please note that only the wavelets of types 1 (orthogonal), 2 (bi-orthogonal) and 3 (discrete Meyer) are possible for MRA. 
Otherwise running errors will occur. 
</p><p>
Since default values are provided with this function, a direct execution of this function without extra input parameters 
will deliver a simple example.
</p><p>
Refer to the description of <a href=\"modelica://Wavelet.Families\">Families</a> for detailed information about the available wavelets.

</p>

</html>"));
    end mra;

    function tuneCoef
      "Tune the wavelet coefficients of all decomposition levels"
      input Real c0[:]= {-0.707, -0.707, -0.707, -0.7069999999999999, -0.707, -0.7069999999999999,
      0.707, 0.0, -1.9993959999999997, -1.9993959999999997, 0.49984899999999993,
      -1.4135729719999994, 3.1805391869999995, 5.6542918879999995, 3.887325672999999}
        "Original wavelet coefficients of all levels";
      input Integer len[:]={13, 7, 4, 2, 2}
        "Lengths of wavelet coefficients for each levels";
      input Real ratio[:]= {1,0.1,2,-1}
        "{rD1, rD2, ..., rDN, rAN}, tuning factors of all levels. ";
      output Real c1[size(c0,1)] "Tuned wavelet coefficients of all levels";

    protected
        Integer k;
        Integer decLevel=size(len,1)-2 "decomposition levels";
        Integer first[decLevel+1]
        "first index of the coefficients for reconstruction";
        Integer last[decLevel+1]
        "last index of the coefficients for reconstruction";

    algorithm
      // initialization
      c1 := c0;

       // coefficient ranges of detail and approximation coefficient vectors
       first := General.cumSumInt(len[2:decLevel+2]) - len[2:decLevel+2] + ones(decLevel+1);
       last := General.cumSumInt(len[2:decLevel+2]);

      for k in 1:size(first,1) loop
        // tuning
        c1[first[k]:last[k]] := c0[first[k]:last[k]] * ratio[k];
      end for;

    annotation (Documentation(info="<html>
<p>
This function supports wavelet MRA. After wavelet (multi-level) decomposition, the coefficients are usually tuned in order to change the signal magnitude at different levels.
This function is able to freely change the magnitude of the coefficients at each level by multiplying them with a given factor, which could be any real value. 
</p><p>
The definition of c[:] and len[:] are the same as in the function <a href=\"modelica://Wavelet.Transform.wavDec\">wavDec()</a>.
</p>

</html>"));
    end tuneCoef;

  annotation (Documentation(info="<html>
<p>
The one dimensional wavelet Multi-Resolution-Analysis (MRA) toolbox is an application function group created on the wavelet 
transform base library. 
Wavelet MRA carries out an N-level wavelet decomposition to obtain the wavelet coefficients in N detail levels and one
approximation level. The coefficients of every single level are then used to reconstruct the signal separately. 
The result is a series of signals, each of which represents a certain frequency range of the original signal. In this way,
the original signal can be observed in different (N+1) resolutions.
</p><p>
This toolbox provides two main functions for MRA. Function <a href=\"modelica://Wavelet.MRA.mra\">mra</a> carries out 
numerical calculation and outputs MRA results in a matrix. Function <a href=\"modelica://Wavelet.MRA.mraGUI\">mraGUI</a> 
provides a graphic user interface to directly display the analysis results. It provides a further possibility to show
the wavelet coefficients in every singles levels.
</p><p>
Only orthogonal, biorthogonal and discrete Meyer wavelets are possible for MRA. Please refer to the description of
<a href=\"modelica://Wavelet.Families\">wavelet families</a> for the detailed information about the available wavelets
in this library.
</p>
</html>"));
  end MRA;

  package Denoising "Application: One dimensional signal denoising"
  function denoisGUI "Graphic user interface for 1D wavelet denoising"
      import Modelica_LinearSystems2.Utilities.Plot.Records.Diagram;
      import Modelica_LinearSystems2.Utilities.Plot.Records.Curve;

    input Real u[:] = {-1.3853,    1.5788,    2.2041,    3.9822,    3.3831,    6.7019,    5.1197,    5.9217,
      4.2061,    3.5389,    3.4177,    2.1313,   -0.7999,   -2.4885,   -3.2673,   -4.6569,
     -4.0413,   -5.1036,   -5.3449,   -6.2969,   -3.7925,   -3.4697,   -3.0802,   -0.8552,
     -1.5134,    1.0541,    2.2248,    3.3398,    3.1799,    4.8989,    4.4210,    4.0804,
      3.4891,    5.1387,    2.3952,    1.9965,   -0.2529,   -0.6175,   -4.8076,   -3.8146,
     -5.1278,   -5.4474,   -5.0879,   -5.3228,   -2.8489,   -4.8475,   -3.7284,   -4.0229,
     -1.7457,   -1.3380} "Signal to be analyzed"                                                                   annotation(Dialog(group="Input data"));
    input Records.wavletDefinition wd
        "Wavelet definition. Note: Only (bi)orthogonal wavelets can be used for MRA. Valid parameters are wavID (<=7), Nd, Nr and range (dMeyer)"
      annotation(Dialog(group="Input data"));
    input Records.denoisParameters dp
        "Data and parameters for wavelet denoising"                               annotation(Dialog(group="Input data"));

    output Real y[:] "Denoised signal" annotation(Dialog(enable = false));
    output Real noise[:] "Noise removed from signal" annotation(Dialog(enable = false));
    output Real c[:]
        "Original wavelet coefficients of all levels, details see wavDec()"
                                                                          annotation(Dialog(enable = false));
    output Real cth[:]
        "Wavelet coefficients of all levels after applying thresholds"                annotation(Dialog(enable = false));
    output Integer len[:]
        "Lengths of coefficients for each levels, details see wavDec()"
                                                                       annotation(Dialog(enable = false));

    protected
    Records.wavFuncOut wOut; // output data of the function wavFunc()
    Integer levels "decomposition levels, <= 12";
    Real thresholds[:]
        "thresholds for all levels with thresholds[1] for the first detail level";
    Integer nSig=size(u,1) "signal length";
    Diagram diagrams[12+2,2]
        "(N+2)x2 martix of diagrams. Note: Non-explicit definition of diagram dimensions does NOT work!";
    String str1 "temp string";
    String str2 "temp string";
    Integer k "temp integer";
    Integer j "temp integer";
    Real ra[:] "temp real array";

  algorithm
    // check the input parameters
    levels := dp.decLevel;
    if levels > 12 then
      levels :=12;
    end if;

    // get the wavelet filters
    wOut:=Families.wavFunc(wd);

    // get the thresholds
    thresholds := {dp.thD1, dp.thD2, dp.thD3, dp.thD4, dp.thD5, dp.thD6, dp.thD7, dp.thD8, dp.thD9, dp.thD10, dp.thD11, dp.thD12};

    // wavelet denoising
    (y, noise, c, cth, len) :=denois(u,levels, wOut.lod, wOut.hid, wOut.lor, wOut.hir, thresholds, dp.sorh, dp.thMethod);

    // --------------------------------------------------------------
    // generate the diagram data
    // --------------------------------------------------------------
    // the original data
    str1 := "Original data and wavelet coefficients";
    diagrams[1,1] := Diagram(curve={Curve(x=(1:nSig), y=u, legend="Original data")},
                            xLabel="Point",yLabel="u",
                            heading=str1);

    // the denoised data
    str1 := "Denoised data and thresholded coefficients";
    diagrams[1,2] := Diagram(curve={Curve(x=(1:nSig), y=y, legend="Denoised data")},
                            xLabel="Point",yLabel="y",
                            heading=str1);

    // The coefficients
    for k in levels+1:-1:1 loop  // loop for all levels
      for j in 1:2 loop // loop for original and denoised data

          if j==1 then // original data
            ra := Transform.wavCoef1(c, len, k);
            str2 :="Original coefficients";
          else // denoised data
            ra := Transform.wavCoef1(cth, len, k);
            str2 :="Thresholded coefficients";
          end if;

        // for the approximation level the label is different as that for the detail levels
        if k==levels+1 then // for approximation
          str1 := "A" + String(levels);
        else  // detail levels
          str1 := "D" + String(k);
        end if;

        // generate the diagram data
        diagrams[levels+3-k,j] := Diagram(curve={Curve(x=(1:size(ra,1)), y=ra, legend=str2)},
                                   xLabel="Point",yLabel=str1);

      end for;
    end for;

    // --------------------------------------------------------------
    // display signals and wavelet coefficients
    Modelica_LinearSystems2.Utilities.Plot.diagramMatrix(diagrams[1:levels+2,:]);

    // display noise data
    Modelica_LinearSystems2.Utilities.Plot.diagram(Diagram(curve={Curve(x=(1:size(noise,1)), y=noise, legend="Noise")},
              xLabel="Point",yLabel="noise",heading="Noise"));
  annotation (Documentation(info="<html>
<p>
This function is a modified version of the function <a href=\"modelica://Wavelet.Families\">denois</a>. 
It provides graphic user interface (GUI) for parameter input and displaying denoising results. 
Instead of filter banks, the wavelet is specified with its name. The original data, the decomposed data in 
different levels, the denoised data in all levels and the denoised data are displayed with diagrams. 
In addition, the noise that is removed from the signal is displayed in a diagram, too.
</p><p>
The default input values provide an example to show the functionality of this function. 
</p>
</html>"));
  end denoisGUI;

    function denois "1D wavelet denoising"

      input Real u[:]={-1.3853,    1.5788,    2.2041,    3.9822,    3.3831,    6.7019,    5.1197,    5.9217,
        4.2061,    3.5389,    3.4177,    2.1313,   -0.7999,   -2.4885,   -3.2673,   -4.6569,
       -4.0413,   -5.1036,   -5.3449,   -6.2969,   -3.7925,   -3.4697,   -3.0802,   -0.8552,
       -1.5134,    1.0541,    2.2248,    3.3398,    3.1799,    4.8989,    4.4210,    4.0804,
        3.4891,    5.1387,    2.3952,    1.9965,   -0.2529,   -0.6175,   -4.8076,   -3.8146,
       -5.1278,   -5.4474,   -5.0879,   -5.3228,   -2.8489,   -4.8475,   -3.7284,   -4.0229,
       -1.7457,   -1.3380} "Data to be denoised";
      input Integer decLevel=2 "Wavelet decomposition levels";
      input Real lod[:]={0.7071067811865476, 0.7071067811865476}
        "Wavelet low pass filter for decomposition";
      input Real hid[:]={-0.7071067811865476, 0.7071067811865476}
        "Wavelet high pass filter for decomposition";
      input Real lor[:]={0.7071067811865476, 0.7071067811865476}
        "Wavelet low pass filter for reconstruction";
      input Real hir[:]={0.7071067811865476, -0.7071067811865476}
        "Wavelet high pass filter for reconstruction";
      input Real th[:]= {-1, -1}
        "Thresholds for all decomposition levels with th[1] for the first detail level (highest frequency components)";
      input Boolean sorh = true
        "Soft (false) or hard (true) thresholding for all decomposition levels";
      input Types.threshMethod thMethod= 1
        "Method for automatic threshold calculation";

      output Real y[:] "Denoised data";
      output Real noise[:] "Noise removed from the data, = u-y";
      output Real c[:] "Wavelet coefficients before denoising";
      output Real cth[:] "Wavelet coefficients after denoising";
      output Integer len[:] "Length of wavelet coefficients of all levels";

    protected
      Integer n = size(u,1) "length of input signal";
      Real autoTh "automatically calculated threshold";
      Integer first[decLevel] "first indices of detail coefficients";
      Integer last[decLevel] "last indices of detail coefficients";

    algorithm
       //wavelet decomposition
       (c, len) := Transform.wavDec(u, decLevel, lod, hid);

       // automatic thresholding value
       if min(th) < 0 then
            autoTh :=thSelect(u, thMethod);
       end if;

       // coefficient ranges for detail levels
       first := General.cumSumInt(len[2:decLevel+1]) - len[2:decLevel+1] + ones(decLevel);
       last := General.cumSumInt(len[2:decLevel+1]);

       //wavelet coefficients thresholding
       cth := c;

       for k in 1:decLevel loop
         if th[k]<0 then
           //thresholding with calculated value
           cth[first[k]:last[k]] := thApply(c[first[k]:last[k]],autoTh,sorh); //apply threshold
         else
           cth[first[k]:last[k]] := thApply(c[first[k]:last[k]],th[k],sorh); //apply threshold
         end if;
       end for;

       //wavelet reconstruction
       y := Transform.wavRec(cth, len, lor,hir);

       noise := u - y;
    annotation (Documentation(info="<html>
<p>
This function carries out the wavelet denoising calculation for a data. The wavelet filter bank must be available for calling
this function. This has to be obtained by wavelet filter functions.
</p><p>
It has to be noted that the thresholds calculated by thSelect in this toolbox 
are valid only if the noise has <b>normal (Gaussian) distribution</b> with <b>zero mean </b> and <b>standard deviation of one</b>!
Otherwise, the data have to be pre-conditioned or the thresholds have to be manually selected by the user. 
</p><p>
The manual selection of the thresholds has the possibility that the thresholds for the levels can be different. The thresholds
are given in the input parameter, th[:], each element of which is used for one level. If a threshold for a level is less than
zero, the default threshold automatically calculated by <a href=\"modelica://Wavelet.Denoising.thSelect\">thSelect</a> will be used.
</p>
</html>"));
    end denois;

    function thSelect "Calculates threshold value for denoising"

      import Modelica.Math.*;

      input Real u[:]={-1.3853,    1.5788,    2.2041,    3.9822,    3.3831,    6.7019,    5.1197,    5.9217,
        4.2061,    3.5389,    3.4177,    2.1313,   -0.7999,   -2.4885,   -3.2673,   -4.6569,
       -4.0413,   -5.1036,   -5.3449,   -6.2969,   -3.7925,   -3.4697,   -3.0802,   -0.8552,
       -1.5134,    1.0541,    2.2248,    3.3398,    3.1799,    4.8989,    4.4210,    4.0804,
        3.4891,    5.1387,    2.3952,    1.9965,   -0.2529,   -0.6175,   -4.8076,   -3.8146,
       -5.1278,   -5.4474,   -5.0879,   -5.3228,   -2.8489,   -4.8475,   -3.7284,   -4.0229,
       -1.7457,   -1.3380} "Signal vector";
      input Types.threshMethod thID= 1 "Threshold calculation method";
      output Real th "Threshold value";

    protected
      Integer n "length of input signal";
      //internal variables for thID=2
      Real sx2[:];
      Real sx2cum[:];
      Real risks[:];
      Integer ind[:];
      Integer best;
      //internal variables for thID=3
      Real hthr;
      Real eta;
      Real crit;

    algorithm
      n := size(u,1);

      if thID == 1 then
        //fixed form thresholding
        th := sqrt(2*log(n));

      elseif thID == 2 then
        //adaptive selection using principle of Stein's unbiased risk estimate
        sx2 := Vectors.sort(abs(u)).^ 2;
        sx2cum := zeros(n);
        sx2cum[1]:= sx2[1];
        for i in 2:n loop
          sx2cum[i] := sx2[i]+sx2cum[i-1];
        end for;
        risks := (n*ones(n)-(2*(1:n)) + (sx2cum + (n-1:-1:0).*sx2))/n;
        (,ind) := Vectors.sort(risks);
        best := ind[1];
        th := sqrt(sx2[best]);

      elseif thID == 3 then
        //heuristic variant of 2
        hthr := sqrt(2*log(n));
        eta := ((Vectors.norm(u))^2-n)/n;
        crit := (log(n)/log(2))^(1.5)/sqrt(n);
        if eta < crit then
          th := hthr;
        else
          th := min(thSelect(u, 2), hthr);
        end if;

      elseif thID == 4 then
        //minimax thresholding
        if n <= 32 then
          th := 0;
        else
          th := 0.3936 + 0.1829*(log(n)/log(2));
        end if;
      end if;
    annotation (Documentation(info="<html>
<p>
This function analyzes the input data and calculates a common threshold for wavelet denoising for all detail levels
using the specified methods. Four methods are available:
</p><p>
1. Fixed form: The threshold is calculated as  sqrt(2*log(size(u,1))).
</p><p>
2. SURE: The threshold is calculated based on Stein's Unbiased Risk Estimate (SURE) of risk, which is a quadratic loss function.
</p><p>
3. Heuristic SURE: This is a heuristic version of the SURE method. It is actually a combination of the first two methods. For the signal with strong
noises, the fixed form method is used. Otherwise, SURE method id used.
</p><p>
4. Minimal maximum: The threshold is estimated to realize the minimum of the maximum mean square error of the observed data. This method is derived 
from statistic theory.
</p><p>
It has to be noted that the thresholds calculated by thSelect in this toolbox 
are valid only if the noise has <b>normal (Gaussian) distribution</b> with <b>zero mean </b> and <b>standard deviation of one</b>!
Otherwise, it is recommended that the data are to be pre-conditioned or the thresholds be manually selected.
</p>
</html>"));
    end thSelect;

    function thApply "Apply soft or hard thresholding"

      input Real u[:]={0,1,2,3,4,5,6,7,8,9} "Input vector";
      input Real th=5 "Threshold value";
      input Boolean sorh = true "Soft (false) or hard (true) thresholding";

      output Real y[:] = zeros(size(u,1)) "Vector after thresholding";

    protected
                Real tmp[:] "temp vector for soft thresholding";
                Real absu[:] = abs(u) "absolute values of input vector";
    algorithm

      if sorh then
        // hard thresholding
        for i in 1:size(u,1) loop
          if absu[i]>th then
            y[i] := u[i];
          end if;
        end for;
      else
        // soft thresholding
        tmp := abs(u) - th*ones(size(u, 1));
        tmp := (tmp + abs(tmp))/2;
        y := sign(u) .* tmp;
      end if;

    end thApply;
  annotation (Documentation(info="<html>
<p>
The wavelet denoising toolbox is an application function group created on the wavelet transform base library. 
Wavelet denoising is based on the consideration that the noise energy in the data is significantly lower than the 
signal energy, or the part of the data that contains useful information. In this case, the noise is mainly represented 
with the wavelet coefficients that have small magnitudes. If these coefficients are removed (set to be zeros) and the signal 
 is reconstructed, the noise will be eliminated from the original data. Therefore, wavelet denoising is only suitable for 
the cases where the signal-to-noise ratio is high.
</p><p>
The thresholds to decide which coefficients are to be removed are usually empirically determined based on the data to 
be analyzed. However, there are methods (the function <a href=\"modelica://Wavelet.Denoising.thSelect\">thSelect</a>) to automatically select the thresholds by analyzing the input data if the noise 
has normal distribution property. It has to be noted that the thresholds calculated by thSelect in this toolbox 
are valid only if the noise has <b>normal (Gaussian) distribution</b> with <b>zero mean </b> and <b>standard deviation of one</b>!
Otherwise, the data have to be pre-conditioned or the thresholds have to be manually selected by the user. 
</p><p>
The same data processing for denoising could also be used for data compression. The idea is based on such consideration: 
The wavelet coefficients with small magnitudes contain little information about the original signal. By removing these 
coefficients, the data amount can be greatly reduced while not much information is lost. If a suitable wavelet is used,
most coefficients might have very small magnitudes. Thus, a large compression ratio can be achieved.
</p><p>
Only orthogonal, biorthogonal or discrete Meyer wavelets can be used for wavelet denoising.
</p>
</html>"));
  end Denoising;

  package Transform "Functions for wavelet transform"
    function cwt
      "One-dimensional continuous wavelet transform with a given wavelet function"

      import Modelica.Math.Vectors.*;
      import Wavelet.General.*;

      input Real[:] u = {1.0000, 1.0000, 0.9999, 0.9996, 0.9987, 0.9968, 0.9933, 0.9877,
     0.9790, 0.9665, 0.9491, 0.9257, 0.8954, 0.8568, 0.8091, 0.7510,
     0.6818, 0.6008, 0.5075, 0.4020, 0.2847, 0.1566, 0.0194,-0.1245,
    -0.2720,-0.4192,-0.5614,-0.6932,-0.8089,-0.9022,-0.9672,-0.9981,
    -0.9901,-0.9397,-0.8452,-0.7074,-0.5297,-0.3186,-0.0837, 0.1623,
     0.4043, 0.6256, 0.8087, 0.9371, 0.9970, 0.9785, 0.8779, 0.6985,
     0.4518, 0.1572,-0.1588,-0.4646,-0.7266,-0.9127,-0.9970,-0.9634,
    -0.8096,-0.5491,-0.2111, 0.1616, 0.5174, 0.8026, 0.9702, 0.9879,
     0.8457, 0.5600, 0.1736,-0.2491,-0.6323,-0.9017,-1.0000,-0.9003,
    -0.6148,-0.1960, 0.2705, 0.6827, 0.9443, 0.9880, 0.7949, 0.4037,
    -0.0939,-0.5730,-0.9055,-0.9959,-0.8104,-0.3923, 0.1445, 0.6439,
     0.9522, 0.9671, 0.6743, 0.1591,-0.4140,-0.8520,-0.9988,-0.7935,
    -0.2997, 0.3091, 0.8072, 1.0000} "Signal to be analyzed";
      input Real[:] scales = {1,4, 16} "Scales to be transformed";
      input Real[:] w = {1,-1} "Wavelet function for the transform";
      input Real[:] x_w = {0,1} "Regular grid of the wavelet function";
      output Real[size(scales,1),size(u,1)] coefs "A two dimensional matrix of wavelet coefficients. 
  Each column corresponds to one element of the input parameter, scales.";

    protected
      Real w1[size(w,1)];
      Integer len_u =  size(u, 1);
      Real step_u =  1;

      Integer n_scales = size(scales, 1);
      Real a;
      Real a_u;

      Real x_w_max;
      Real[size(x_w,1)] x_w_new "aligned to the origin";
      Real step_w;

      Integer j[:];
      Integer n_j;
      Real f[:];

    algorithm
        x_w_new := x_w .- x_w[1];
        step_w :=x_w_new[2] - x_w_new[1];

        w1 := step_w .* General.cumSum(w);  // integration

        x_w_max := x_w_new[size(x_w, 1)];

        for i in 1:n_scales loop
            a := scales[i];
            a_u := a/step_u;
            j := integer((0:a_u*x_w_max) / (a_u*step_w));

            n_j := size(j,1);
            f := fill(0, n_j);

            for k in 1:n_j loop
                f[n_j+1-k] := w1[(j[k]+1)];
            end for;

            coefs[i,:] := -sqrt(a).*midVector(diff(wavConv(u, f)), len_u);

        end for;

    annotation (Documentation(info="<html>
<p>
One-dimensional CWT transforms a vector to a two dimensional image. The first dimension (horizontal axis) 
is the same as the original data. The second dimension (vertical axis) is the scale, 
similar to frequency. The scales are defined by the input parameter, scales. Larger scale values correspond to lower frequencies.
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example with a simple saw-tooth signal.</p>

</html>"));
    end cwt;

    function cwtn
      "One-dimensional continuous wavelet transform with a given wavelet name"
        input Real[:] u = {1.0000, 1.0000, 0.9999, 0.9996, 0.9987, 0.9968, 0.9933, 0.9877,
     0.9790, 0.9665, 0.9491, 0.9257, 0.8954, 0.8568, 0.8091, 0.7510,
     0.6818, 0.6008, 0.5075, 0.4020, 0.2847, 0.1566, 0.0194,-0.1245,
    -0.2720,-0.4192,-0.5614,-0.6932,-0.8089,-0.9022,-0.9672,-0.9981,
    -0.9901,-0.9397,-0.8452,-0.7074,-0.5297,-0.3186,-0.0837, 0.1623,
     0.4043, 0.6256, 0.8087, 0.9371, 0.9970, 0.9785, 0.8779, 0.6985,
     0.4518, 0.1572,-0.1588,-0.4646,-0.7266,-0.9127,-0.9970,-0.9634,
    -0.8096,-0.5491,-0.2111, 0.1616, 0.5174, 0.8026, 0.9702, 0.9879,
     0.8457, 0.5600, 0.1736,-0.2491,-0.6323,-0.9017,-1.0000,-0.9003,
    -0.6148,-0.1960, 0.2705, 0.6827, 0.9443, 0.9880, 0.7949, 0.4037,
    -0.0939,-0.5730,-0.9055,-0.9959,-0.8104,-0.3923, 0.1445, 0.6439,
     0.9522, 0.9671, 0.6743, 0.1591,-0.4140,-0.8520,-0.9988,-0.7935,
    -0.2997, 0.3091, 0.8072, 1.0000} "Signal to be analyzed"
                               annotation(Dialog(group="Input data"));
        input Real[:] scales = {1,4,16} "Scales for the transformation"
                                                                       annotation(Dialog(group="Input data"));
        input Records.wavletDefinition wIn "Wavelet definition"
                            annotation(Dialog(group="Input data"));

        output Real[size(scales,1),size(u, 1)] coefs "A two dimensional matrix of wavelet coefficients. 
  Each column corresponds to one element of the input parameter, scales.";

    protected
        Records.wavFuncOut wOut "output data of the function wavFunc()";

        Real step "step of the wavelet";
    algorithm

        // get the wavelet
        wOut := Families.wavFunc(wIn);

        // cwt
        coefs := Transform.cwt(u, scales, wOut.psi1, wOut.x);
    annotation (Documentation(info="<html>
<p>
This function carries out one-dimensional CWT transform, same as the function, cwt(). 
The only difference is that cwtn() receives the wavelet name as input. 
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example with a simple saw-tooth signal.</p>

</html>"));
    end cwtn;

    function dwt
      "One-dimensional discrete wavelet transform (one-level decomposition)"
          input Real u[:]={1,2,1,2,1,2,3,4,1} "Signal to be analyzed";
          input Real lod[:]={0.707, 0.707} "Low pass filter for decomposition";
          input Real hid[size(lod,1)]={-0.707, 0.707}
        "High pass filter for decomposition";
          output Real ca[:]
        "Approximation coefficients, size(ca,1) = ceil(size(u,1)/2)";
          output Real cd[size(ca,1)] "Detail coefficients";

    protected
          Real y[:] "convolution result";
          Integer nu "length of the input signal";
          Integer nf "length of the filters";

    algorithm
          // the vector lengthes
          nu :=size(u, 1);
          nf :=size(lod, 1);

          // transformation
          y :=Wavelet.General.wavConv(u, lod);      // approximation coefficients
          ca :=y[2:2:size(y,1)];  // down-sampling

          y :=Wavelet.General.wavConv(u, hid);      // detail coefficients
          cd :=y[2:2:size(y,1)];  // down-sampling

          annotation (Inline=true, Documentation(info="<html>
<p>
This function calculates discrete wavelet transform (one-level wavelet decomposition) of the input signal, u. 
Before transformation, the original signal is extended at both ends with zero padding. 
The length of the two result vectors is floor((size(u,1)+size(lod,1)-1)/2).
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example with a simple saw-tooth signal.
</p>
</html>"));
    end dwt;

    function idwt
      "One-dimensional inverse discrete wavelet transform (one-level reconstruction)"
          input Real ca[:]={2.121, 2.121, 2.121, 4.95, 0.707}
        "Approximation coefficients";
          input Real cd[size(ca,1)]={-0.707, -0.707, -0.707, -0.707, 0.707}
        "Detail coefficients";
          input Real lor[:]={0.707, 0.707} "Low pass filter for reconstruction";
          input Real hir[size(lor,1)]={0.707, -0.707}
        "High pass filter for reconstruction";
          input Integer ny=0
        "Length of the result, y (the central part will be extracted). If ny == 0, size(y,1) = 2*size(ca,1) - size(lor,1) +2";
          output Real y[:] "Reconstructed signal";

    protected
          Integer nc = size(ca,1) "input data length";
          Integer nf = size(lor, 1) "filter length";

    algorithm
          // inverse transformation
          y:=General.wavConv(Wavelet.General.upsample(ca), lor) + General.wavConv(Wavelet.General.upsample(cd), hir);
          if ny == 0 then
            y := General.midVector(y, 2*nc - nf +2);
          else
            y := General.midVector(y, ny);
          end if;

          annotation (Inline=true, Documentation(info="<html>
<p>
This function calculates inverse discrete wavelet transform (one-level wavelet reconstruction) using the input approximation 
and detail wavelet coefficients, ca and cd. 
The coefficient vectors are extended at both ends with zero padding before calculation. 
The length of the reconstructed data is usually (but not always) equal to the length of the original signal if the length is not specified.
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example to generate a simple saw-tooth signal.
</p>
</html>"));
    end idwt;

    function wavDec "Wavelet multilevel decomposition"
          input Real u[:]={1,2,1,2,1,2,3,4,1,2,3,4,1} "Signal to be analyzed";
          input Integer N=3 "Decomposition level";
          input Real lod[:]={0.707, 0.707} "Low pass filter for decomposition";
          input Real hid[size(lod,1)]={-0.707, 0.707}
        "High pass filter for decomposition";
          output Real c[:]
        "Detail and approximation coefficients. For data structure see the information section";
          output Integer len[N+2]
        "Lengths of the coefficient vectors of all levels. For data structure see the information section";

    protected
          Real ca[:] "temp variable for approximation coefficients";
          Real cd[:] "temp variable for detail coefficients";
          Integer k "level counter";
          Integer i "debug";

    algorithm
          // first level dwt
          (ca,cd) := dwt(u, lod, hid);
          c := cd;
          len[1] := size(u,1);
          len[2] := size(cd,1);

          // higher level dwt
          for k in 2:N loop
            (ca,cd) := dwt(ca, lod, hid);
            c := cat(1, c, cd);  // attach the new detail coef at the beginning of c
            len[k+1] := size(cd,1);

          end for;

          // approx coef of highest level (lowest frequency)
          c := cat(1, c, ca);
          len[N+2] := size(ca,1);
    annotation (Documentation(info="<html>
<p>
This function carries out multi-level wavelet decomposition. </p><p>
Parameter, len, is a vector with (N+2) elements, where N is the number of levels. 
It contents the length information of the coefficients in all levels. 
len[1] is the length of the original signal vector; 
len[2] is the length of the detail coefficient vector of the first decomposition level; 
len[3] is the length of the detail coefficient vector of the second decomposition level; and so forth.
Finally, len[N+2] contains the length of the approximation coefficient vector.
Parameter, c, contains the wavelet coefficients of all levels, firstly the first detail level, 
then the second detail level and so on, and finally the approximation coefficients.
This is described again in detail as follows:
</p>
<h4><font color=\"#008000\">Definition of the outputs</font></h4>

<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>Data</b></td>      <td><b>Length</b></td>      <td><b>Location</b></td>      </tr>
        
  <tr><td>Original signal</td>                    <td>len[1]</td>   <td>u</td></tr>
  <tr><td>First level detail coefficients</td>   <td>len[2]</td>   <td>c[1 : len[2]]</td></tr>
  <tr><td>Second level detail coefficients</td>   <td>len[3]</td>   <td>c[len[2]+1 : len[2]+len[3]]</td></tr>
  <tr><td>...</td>                                <td>...</td>       <td>...</td></tr>
  <tr><td>k-th level detail coefficients</td>   <td>len[k+1]</td>   <td>c[sum(len[2:k])+1 : sum(len[2:k+1])]</td></tr>
  <tr><td>...</td>                              <td>...</td>        <td>...</td></tr>
  <tr><td>N-th level detail coefficients</td>   <td>len[N+1]</td>   <td>c[sum(len[2:N])+1 : sum(len[2:N+1])]</td></tr>
  <tr><td>Approximation coefficients</td>        <td>len[N+2]</td>   <td>c[sum(len[2:N+1])+1 : sum(len[2:N+2])]</td></tr>

</table>
<p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example with a simple saw-tooth signal.
</p>
</html>"));
    end wavDec;

    function wavRec "Wavelet multilevel reconstruction"
        input Real c[:]= {-0.707, -0.707, -0.707, -0.7069999999999999, -0.707, -0.7069999999999999,
      0.707, 0.0, -1.9993959999999997, -1.9993959999999997, 0.49984899999999993,
      -1.4135729719999994, 3.1805391869999995, 5.6542918879999995, 3.887325672999999}
        "Detail and approximation coefficients. For data structure see the information section";
        input Integer len[:]={13, 7, 4, 2, 2}
        "Coefficient vector lengths of all levels. For data structure see the information section";
        input Real lor[:]={0.707, 0.707} "Low pass filter for reconstruction";
        input Real hir[size(lor,1)]={0.707, -0.707}
        "High pass filter for reconstruction";
        output Real u[:] "Reconstructed signal";

    protected
        Real cd[:] "temp variable for detail coefficients";
        Integer decLevel= size(len,1)-2 "number of levels";
        Integer first[decLevel+1]
        "first index of the coefficients for reconstruction";
        Integer last[decLevel+1]
        "last index of the coefficients for reconstruction";
        Integer k "level counter";

    algorithm
    /* debug: display a vector in error message window
    Modelica.Utilities.Streams.print("size(c,1) = " + String(size(c,1)));
    for i in 1:size(c,1) loop
      Modelica.Utilities.Streams.print(String(c[i]));
    end for;
// debug end */

       // coefficient ranges of detail and approximation coefficient vectors
       first := General.cumSumInt(len[2:decLevel+2]) - len[2:decLevel+2] + ones(decLevel+1);
       last := General.cumSumInt(len[2:decLevel+2]);

      // get the approximation coefficients
      u := c[first[decLevel+1] : last[decLevel+1]];

      // reconstruction from higher to lower levels
        for k in decLevel:-1:1 loop
            cd := c[first[k] : last[k]];
            u := idwt(u, cd, lor, hir, len[k]);
        end for;

    annotation (Documentation(info="<html>
<p>
This function transforms the wavelet coefficients, which are obtained with multi-level wavelet decomposition, back to 
the original signal. The definitions of the input parameters, c and len, are same as the function <a href=\"modelica://Wavelet.Transform.wavDec\">wavDec()</a>.
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example to generate a simple saw-tooth signal.
</p>

</html>"));
    end wavRec;

    function wavRec1
      "Wavelet reconstruction using coefficients from a single level"
      input Real c[:]={-0.707, -0.707, -0.707, -0.7069999999999999, -0.707, -0.7069999999999999,
      0.707, 0.0, -1.9993959999999997, -1.9993959999999997, 0.49984899999999993,
      -1.4135729719999994, 3.1805391869999995, 5.6542918879999995, 3.887325672999999}
        "Wavelet coefficients of all levels";
      input Integer len[:]={13, 7, 4, 2, 2}
        "Lengths of coefficients of all levels";
      input Real lor[:]={0.707, 0.707} "Low pass filter for reconstruction";
      input Real hir[size(lor,1)]={0.707, -0.707}
        "High pass filter for reconstruction";
      input Integer level= 2 "The level to be reconstructed. 1..N are for the detail levels, 
    where N is the maximum decomposition level. 1 stands for the most detail level (highest frequency).
     N+1 is for the approximation level.";
      output Real y[:] "Reconstructed signal";

    protected
      Real c1[size(c,1)] "modified coefficients for reconstruction";
      Integer first=1 "first index of the coefficients for reconstruction";
      Integer last=1 "last index of the coefficients for reconstruction";

    algorithm
      // check the input
      if level < 1 or level > size(len,1)-1 then
        return;
      end if;

      // coefficient range
       first := sum(len[2:level+1]) - len[level+1] + 1;
       last  := sum(len[2:level+1]);

      // coefficients for reconstruction
      c1 := zeros(size(c,1));
      c1[first:last] :=  c[first:last];

      // reconstruction
      y := Transform.wavRec(c1, len, lor, hir);

    annotation (Documentation(info="<html>
<p>
This function carries out a modified wavelet reconstruction: 
It reconstructs the signal with the wavelet coefficients from only one single level. 
All other coefficients are set zero.
This operation is useful for some wavelet analysis, e.g. multi-resolution analysis and wavelet de-noising.
</p><p>
The definitions of input parameters, c and len, are same as the function <a href=\"modelica://Wavelet.Transform.wavDec\">wavDec()</a>.
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example.
</p>

</html>"));
    end wavRec1;

  function wavCoef1 "Extract the wavelet coefficients of a single level"
    input Real c[:]={-0.707, -0.707, -0.707, -0.7069999999999999, -0.707, -0.7069999999999999,
    0.707, 0.0, -1.9993959999999997, -1.9993959999999997, 0.49984899999999993,
    -1.4135729719999994, 3.1805391869999995, 5.6542918879999995, 3.887325672999999}
        "Wavelet coefficients of all levels";
    input Integer len[:]={13, 7, 4, 2, 2}
        "Lengths of coefficient vectors for all levels";
    input Integer level = 2 "The level to be reconstructed. 1..N are for the detail levels, 
    where N is the maximum decomposition level. 1 stands for the most detail level (highest frequency).
     N+1 is for the approximation level.";
    output Real c1[:] "Extracted coefficients";

    protected
    Integer first=1 "first index of the coefficients for reconstruction";
    Integer last=1 "last index of the coefficients for reconstruction";

  algorithm
    // check the input
    if level < 1 or level > size(len,1)-1 then
      return;
    end if;

    // coefficient range
    first := sum(len[2:level+1]) - len[level+1] + 1;
    last  := sum(len[2:level+1]);

    // coefficients for reconstruction
    c1 :=  c[first:last];

  annotation (Documentation(info="<html>
<p>
This function extracts the coefficients of a specified decomposition level from the data structure of the wavelet coefficients. 
This operation is useful if only the coefficients of a specific level are interested.
</p><p>
The definitions of the input parameters, c and len, are same as the function <a href=\"modelica://Wavelet.Transform.wavDec\">wavDec()</a>.
</p><p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for the detailed information about the available wavelets.
</p><p>
The default input values provide a quick example.
</p>

</html>"));
  end wavCoef1;

  annotation (Documentation(info="<html>
<p>
This library includes the following one-dimensional wavelet transforms and inverse transforms:
</p><p>
1. Continuous wavelet transform: cwt() and cwtn();</p><p>
2. Discrete wavelet transform: dwt();</p><p>
3. Inverse discrete wavelet transform: idwt();</p><p>
4. Multi-level wavelet decomposition: wavDec();</p><p>
5. Multi-level wavelet reconstruction: wavRec();</p><p>
6. Single level reconstruction based on multi-level wavelet coefficients: wavRec1();
</p></html>"));

  end Transform;

  package Families "Functions about wavelet families"
  function wavFunc
      "Returns the filter bank and functions of a specific wavelet family"
      import Modelica_LinearSystems2.Math.Complex;

      //input Records.wavFuncIn wIn
      input Records.wavletDefinition wIn
        "Input parameters for wavelet filters and functions";
      output Records.wavFuncOut wOut "Output data of the function wavFunc()";

    protected
    Real lod1[:];  // wavelet filters for biorthogonal wavelets
    Real hid1[:];
    Real lor1[:];
    Real hir1[:];
    Real lod2[:];
    Real hid2[:];
    Real lor2[:];
    Real hir2[:];

    Real mul;
    Real es=wIn.range; // effective support
    Complex cpsi[:];

  algorithm
    // Haar wavelet -----------------------
    if wIn.wavID==1 then
      wOut.wavName :="Haar";
      wOut.wavType :=1;
      (wOut.F1, wOut.lod, wOut.hid, wOut.lor, wOut.hir):= Families.wavHaar();
      (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, wOut.lor, wOut.hir, wIn.refinement);  // scaling and wavelet functions

    // Daubechies wavelet -----------------------
    elseif wIn.wavID==2 then
      wOut.wavName := "Daubechies";
      wOut.wavType :=1;
      (wOut.F1, wOut.lod, wOut.hid, wOut.lor, wOut.hir):= Families.wavDaubechies(wIn.Nd);
      (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, wOut.lor, wOut.hir, wIn.refinement);  // scaling and wavelet functions

    // Symlets wavelet -----------------------
    elseif wIn.wavID==3 then
      wOut.wavName :="Symlets";
      wOut.wavType :=1;
      (wOut.F1, wOut.lod, wOut.hid, wOut.lor, wOut.hir):= Families.wavSymlets(wIn.Nd);
      (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, wOut.lor, wOut.hir, wIn.refinement);  // scaling and wavelet functions

    // Coiflets wavelet -----------------------
    elseif wIn.wavID==4 then
      wOut.wavName :="Coiflets";
      wOut.wavType :=1;
      (wOut.F1, wOut.lod, wOut.hid, wOut.lor, wOut.hir):= Families.wavCoiflets(wIn.Nd);
      (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, wOut.lor, wOut.hir, wIn.refinement);  // scaling and wavelet functions

    // Biorthogonal spline wavelet -----------------------
    elseif wIn.wavID==5 then
      wOut.wavName :="Biorthogonal spline";
      wOut.wavType :=2;
      (wOut.F1, lod1, hid1, lor1, hir1, wOut.F2, lod2, hid2, lor2, hir2):= Families.wavBiorSpline(wIn.Nd, wIn.Nr);

      // filters for transforms
      wOut.lod := lod2;
      wOut.hid := hid1;
      wOut.lor := lor1;
      wOut.hir := hir2;

      // phi and psi functions:
      (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, lor1, hir2, wIn.refinement);  // scaling and wavelet functions for decomposition
      (wOut.x, wOut.phi2, wOut.psi2) := scalingWaveFunc(wIn.wavID, lor2, hir1, wIn.refinement);  // scaling and wavelet functions for reconstruction

      // reverse the psi function form some orders
      mul :=-1;
      if rem(wIn.Nr, 4) == 1 then
        mul :=1;
      end if;
      wOut.psi1 := mul*wOut.psi1;
      wOut.psi2 := mul*wOut.psi2;

    // reverse Biorthogonal spline wavelet -----------------------
    elseif wIn.wavID==6 then
        wOut.wavName :="Reverse Biorthogonal spline";
        wOut.wavType :=2;
        (wOut.F1, lod1, hid1, lor1, hir1, wOut.F2, lod2, hid2, lor2, hir2):= Families.wavRevBiorSpline(wIn.Nd, wIn.Nr);

        // filters for transforms
        wOut.lod := lod2;
        wOut.hid := hid1;
        wOut.lor := lor1;
        wOut.hir := hir2;

        // phi and psi functions:
        (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, lor1, hir2, wIn.refinement);  // scaling and wavelet functions for decomposition
        (wOut.x, wOut.phi2, wOut.psi2) := scalingWaveFunc(wIn.wavID, lor2, hir1, wIn.refinement);  // scaling and wavelet functions for reconstruction

        // reverse the psi function form some orders
        mul :=-1;
        if rem(wIn.Nd, 4) == 1 then
          mul :=1;
        end if;
        wOut.psi1 := mul*wOut.psi1;
        wOut.psi2 := mul*wOut.psi2;

    // discrete Meyer wavelet -----------------------
    elseif wIn.wavID==7 then
        wOut.wavName :="Discrete Meyer";
        wOut.wavType :=3;
        (wOut.F1, wOut.lod, wOut.hid, wOut.lor, wOut.hir):= Families.wavDMeyer();
        (wOut.x, wOut.phi1, wOut.psi1) := scalingWaveFunc(wIn.wavID, wOut.lor, wOut.hir, wIn.refinement);  // scaling and wavelet functions

    // Meyer wavelet -----------------------
    elseif wIn.wavID==8 then
        wOut.wavName :="Meyer";
        wOut.wavType :=3;
        if es <=0 then
          es := 8;  // effective support as default value
        end if;
        (wOut.x, wOut.phi1, wOut.psi1) := Families.wavMeyer(-es, es, integer(2^wIn.refinement));

    // Gaussian wavelet -----------------------
    elseif wIn.wavID==9 then
        wOut.wavName :="Gaussian";
        wOut.wavType :=4;
        if es <=0 then
          es := 5;  // effective support as default value
        end if;
        (wOut.x, wOut.psi1) := Families.wavGaussian(-es, es, wIn.Nd, integer(2^wIn.refinement));

    // Mexican hat wavelet -----------------------
    elseif wIn.wavID==10 then
        wOut.wavName :="Mexican hat";
        wOut.wavType :=4;
        if es <=0 then
          es := 5;  // effective support as default value
        end if;
        (wOut.x, wOut.psi1) := Families.wavMexHat(-es, es, integer(2^wIn.refinement));

    // Morlet wavelet -----------------------
    elseif wIn.wavID==11 then
        wOut.wavName :="Morlet";
        wOut.wavType :=4;
        if es <=0 then
          es := 4;  // effective support as default value
        end if;
        (wOut.x, wOut.psi1) := Families.wavMorlet(-es, es, integer(2^wIn.refinement));

    // complex Gaussian wavelet -----------------------
    elseif wIn.wavID==12 then
        wOut.wavName :="Complex Gaussian";
        wOut.wavType :=5;
        if es <=0 then
          es := 5;  // effective support as default value
        end if;
        (wOut.x, cpsi) := Families.wavXGaussian(-es, es, wIn.Nd, integer(2^wIn.refinement));
        wOut.psi1 := Complex.real(cpsi);
        wOut.psi2 := Complex.imag(cpsi);

    // complex Morlet wavelet -----------------------
    elseif wIn.wavID==13 then
        wOut.wavName :="Complex Morlet";
        wOut.wavType :=5;
        if es <=0 then
          es := 4;  // effective support as default value
        end if;
        (wOut.x, cpsi) := Families.wavXMorlet(-es, es, integer(2^wIn.refinement), wIn.fb, wIn.fc);
        wOut.psi1 := Complex.real(cpsi);
        wOut.psi2 := Complex.imag(cpsi);

    // complex Shannon wavelet -----------------------
    elseif wIn.wavID==14 then
        wOut.wavName :="Complex Shannon";
        wOut.wavType :=5;
        if es <=0 then
          es := 20;  // effective support as default value
        end if;
        (wOut.x, cpsi) := Families.wavXShannon(-es, es, integer(2^wIn.refinement), wIn.fb, wIn.fc);
        wOut.psi1 := Complex.real(cpsi);
        wOut.psi2 := Complex.imag(cpsi);

    // complex frequency B-Spline wavelet -----------------------
    elseif wIn.wavID==15 then
        wOut.wavName :="Complex frequency B-Spline";
        wOut.wavType :=5;
        if es <=0 then
          es := 20;  // effective support as default value
        end if;
        (wOut.x, cpsi) := Families.wavXFreqBSpline(-es, es, wIn.Nd, integer(2^wIn.refinement), wIn.fb, wIn.fc);
        wOut.psi1 := Complex.real(cpsi);
        wOut.psi2 := Complex.imag(cpsi);

    end if;

  annotation (Documentation(info="<html>
<p>
This function provides a common entry to access all available wavelets in this library.
Since different wavelet families require different parameters and return different results,
suitable input parameters must be provided. Not applicable parameters for a specific wavelet
are omitted by the function.
</p><p>
Refer to <a href=\"modelica://Wavelet.Families\">Wavelet Families</a> 
for detailed information about the available wavelets.
</p>
</html>"));
  end wavFunc;

  function scalingWaveFunc
      "Generate the scaling and wavelet functions using wavelet filters"
    input Types.wavletID wavID=1
        "Id-Number of the wavelet. Necessary only for discrete Meyer wavelet";
    input Real lof[:]={0.707, 0.707} "Wavelet low pass filter";
    input Real hif[:]={-0.707, 0.707} "Wavelet high pass filter";
    input Integer iteration=6
        "Iteration of the calculations for generating the functions";
    output Real x[:] "Equidistant grid of the functions";
    output Real phi[:] "The scaling function";
    output Real psi[:] "The wavelet function";

    protected
    Real mag = sqrt(2)^iteration "magnitude of the psi function";
    Integer nf = size(lof,1) "filter length";
    Integer nFunc "length of the functions";

  algorithm
    if iteration>=0 then  // if 0, use 1
      // the first iteration
      phi := lof;
      psi := hif;
      nFunc:= nf;

      // the rest iterations
      for k in 2:iteration loop
        // scaling function
        nFunc:= 2*size(phi,1)+nf-2;
        phi := General.wavConv(General.upsample(phi),lof);
        phi := General.midVector(phi, nFunc);

        // wavelet function
        psi := General.wavConv(General.upsample(psi),lof);
        psi := General.midVector(psi, nFunc);
      end for;

      // correct the magnitude
      phi := mag * phi;
      psi := mag * psi;

      // correct the length. total length = (nf-1)*2^interation+1
      phi := cat(1,{0}, phi, fill(0, integer((nf-1)*2^iteration - nFunc)));
      psi := cat(1,{0}, psi, fill(0, integer((nf-1)*2^iteration - nFunc)));

      // x-axis
      nFunc :=size(phi, 1);
      x := (0:(nFunc-1)) / 2^iteration;

      // correct the function signs acording to the wavelets
      if wavID == Types.wavletID.Symlets or wavID == Types.wavletID.Coiflets or wavID == Types.wavletID.DiscreteMeyer then
          psi := -psi;
      end if;

    end if;

    annotation (Inline=true, Documentation(info="<html>
<p>
The total length of the functions is (size(lof,1)-1)*2^iteration+1. 
The function support (non-zero part) is (2^iteration -1)*nf-2^iteration+2, shorter than total length. 
One 'zero' is then added at the beginning and multiple 'zeros' are added at the end.
</p>
</html>"));
  end scalingWaveFunc;

    function wavHaar "Haar wavelet filters"
        output Real F[:] "Scaling filter";
        output Real lod[:] "High pass filter for decomposition";
        output Real hid[:] "Low pass filter for decomposition";
        output Real lor[:] "High pass filter for reconstruction";
        output Real hir[:] "Low pass filter for reconstruction";

    algorithm
        F :={0.50000000000000,0.50000000000000};

        //----------------------------------------------------
        // filter bank
        (lod, hid, lor, hir) := General.filterBank(F);
    annotation (Documentation(info="<html>
<p>
This function generate the Haar wavelet filter bank. No input parameters are required.
</p>
</html>"));
    end wavHaar;

    function wavDaubechies "Daubechies wavelet filters"
      input Integer order=2 "Order of the wavelet. order<=10";
      output Real F[:] "Scaling filter";
      output Real lod[:] "High pass filter for decomposition";
      output Real hid[:] "Low pass filter for decomposition";
      output Real lor[:] "High pass filter for reconstruction";
      output Real hir[:] "Low pass filter for reconstruction";

    algorithm
      if order==1 then
        F := {5.000000000000000e-01, 5.000000000000000e-01};

      elseif order==2 then
        F := {3.415063509462200e-01, 5.915063509458700e-01, 1.584936490537800e-01,-9.150635094586999e-02};

      elseif order==3 then
        F := {2.352336038927000e-01, 5.705584579173100e-01, 3.251825002637100e-01,-9.546720778426000e-02,-6.041610415535000e-02,
              2.490874986589000e-02};

      elseif order==4 then
        F := {1.629017140256200e-01, 5.054728575456500e-01, 4.461000691231900e-01,-1.978751311791000e-02,-1.322535836843700e-01,
              2.180815023739000e-02, 2.325180053556000e-02,-7.493494665130000e-03};

      elseif order==5 then
        F := {1.132094912917300e-01, 4.269717713527100e-01, 5.121634721301600e-01, 9.788348067375000e-02,-1.713283576913300e-01,
              -2.280056594205000e-02, 5.485132932108000e-02,-4.413400054330000e-03,-8.895935050930000e-03, 2.358713969200000e-03};

      elseif order==6 then
        F := {7.887121600143000e-02, 3.497519070375700e-01, 5.311318799412100e-01, 2.229156614650500e-01,-1.599932994458700e-01,
              -9.175903203003000e-02, 6.894404648720000e-02, 1.946160485396000e-02,-2.233187416548000e-02, 3.916255760300000e-04,
              3.378031181510000e-03,-7.617669025800000e-04};

      elseif order==7 then
        F := {5.504971537285000e-02, 2.803956418130400e-01, 5.155742458183300e-01, 3.321862411056600e-01,-1.017569112317300e-01,
              -1.584175056405400e-01, 5.042323250485000e-02, 5.700172257986000e-02,-2.689122629486000e-02,-1.171997078235000e-02,
              8.874896189620000e-03, 3.037574977600000e-04,-1.273952359060000e-03, 2.501134265800000e-04};

      elseif order==8 then
        F := {3.847781105406000e-02, 2.212336235762400e-01, 4.777430752143800e-01, 4.139082662116600e-01,-1.119286766665000e-02,
              -2.008293163911100e-01, 3.340970462800000e-04, 9.103817842345000e-02,-1.228195052300000e-02,-3.117510332533000e-02,
              9.886079648080000e-03, 6.184422409540000e-03,-3.443859628130000e-03,-2.770022742100000e-04, 4.776148553300000e-04,
              -8.306863060000001e-05};

      elseif order==9 then
        F := {2.692517479416000e-02, 1.724171519247100e-01, 4.276745321702800e-01, 4.647728571727800e-01, 9.418477475112000e-02,
              -2.073758808962800e-01,-6.847677451089999e-02, 1.050341711371400e-01, 2.172633772990000e-02,-4.782363205882000e-02,
              1.774464067300000e-04, 1.581208292614000e-02,-3.339810113240000e-03,-3.027480287150000e-03, 1.306483640180000e-03,
              1.629073360100000e-04,-1.781648795500000e-04, 2.782275679000000e-05};

      elseif order==10 then
        F := {1.885857879640000e-02, 1.330610913968700e-01, 3.727875357426600e-01, 4.868140553661000e-01, 1.988188708844000e-01,
              -1.766681008964700e-01,-1.385549393599300e-01, 9.006372426666000e-02, 6.580149355070000e-02,-5.048328559801000e-02,
              -2.082962404385000e-02, 2.348490704841000e-02, 2.550218483930000e-03,-7.589501167680000e-03, 9.866626824400000e-04,
              1.408843294960000e-03,-4.849739199600000e-04,-8.235450295000000e-05, 6.617718320000000e-05,-9.379207890000000e-06};

      elseif order==11 then
        F := {1.321886471665690e-02, 1.018707676009539e-01, 3.181271742303942e-01, 4.848537683131867e-01, 2.913027988903065e-01,
              -1.147459261776494e-01,-1.939104913955044e-01, 4.669986906776433e-02, 1.059330899181783e-01,-3.286629145225681e-02,
              -4.697931589875463e-02, 2.215725478297892e-02, 1.473674479914390e-02,-1.086456905449246e-02,-2.362343964095948e-03,
              3.484917545118838e-03,-2.182081030939597e-04,-6.314627963033845e-04, 1.761774389539444e-04, 3.849423888144486e-05,
              -2.449063218490632e-05, 3.177931817946260e-06};

      elseif order==12 then
        F := {9.271766518224606e-03, 7.747505450119822e-02, 2.668303750255462e-01, 4.647096733229814e-01, 3.647868272187679e-01,
              -3.165284709804767e-02,-2.235719287136848e-01,-1.681447405735251e-02, 1.290318596728583e-01, 3.789788060971172e-03,
              -6.818780604444188e-02, 7.671493573866041e-03, 2.937765454983170e-02,-8.639889614159433e-03,-9.079834573750825e-03,
              4.745746461045812e-03, 1.590005428332839e-03,-1.541141788352374e-03, 4.628104542809228e-06, 2.748192162496153e-04,
              -6.258185578399792e-05,-1.714136139124106e-05, 9.034669557220284e-06,-1.081217009051121e-06};

      elseif order==13 then
        F := {6.506891026784466e-03, 5.859174744008193e-02, 2.206147151049089e-01, 4.320817360380974e-01, 4.164078086219489e-01,
              6.150819684806527e-02,-2.227194789327674e-01,-8.808905109195202e-02, 1.269087528252618e-01, 5.158268566903442e-02,
              -7.481728432188596e-02,-1.873013184154149e-02, 3.969660494987692e-02, 1.682894519880251e-03,-1.685135918958014e-02,
              2.774645607423246e-03, 5.130476467388951e-03,-1.952966163061167e-03,-9.303219449292584e-04, 6.592541294136748e-04,
              3.482608740058732e-05,-1.167638275787789e-04, 2.169300195922749e-05, 7.383559915721792e-06,-3.323696366957070e-06,
              3.691122216149230e-07};

      elseif order==14 then
        F := {4.568725425913564e-03, 4.409854388945876e-02, 1.802063525433091e-01, 3.919532612956728e-01, 4.463172083044448e-01,
              1.546235261609992e-01,-1.921128176869846e-01,-1.541729875842131e-01, 9.786019420757067e-02, 9.898718291852615e-02,
              -6.134039007693343e-02,-5.059275162364903e-02, 3.905854655118479e-02, 1.907873678051358e-02,-2.134426676668772e-02,
              -3.970439599608046e-03, 9.043537416561069e-03,-5.276565075623354e-04,-2.722105748697295e-03, 7.507289661574468e-04,
              5.006465593832778e-04,-2.735314931244920e-04,-2.954097378351738e-05, 4.861715681157672e-05,-7.309510712951333e-06,
              -3.103990103456402e-06, 1.219755432462842e-06,-1.263698790522041e-07};

      elseif order==15 then
        F := {3.209230545041146e-03, 3.305257150435777e-02, 1.456808713114689e-01, 3.483432664027908e-01, 4.566588509261426e-01,
              2.397109916594859e-01,-1.366159572709547e-01,-2.042708429993033e-01, 4.616201865523868e-02, 1.344540308947747e-01,
              -2.804822242634575e-02,-7.857436740372240e-02, 2.395475819555182e-02, 3.873569879543842e-02,-1.822002561283388e-02,
              -1.471492759183369e-02, 1.066594072433904e-02, 3.606951945677626e-03,-4.587521102137361e-03,-1.709476540125823e-04,
              1.374137564571244e-03,-2.640919052640912e-04,-2.542510225676878e-04, 1.102354651860462e-04, 1.823819247783242e-05,
              -1.989324456685357e-05, 2.377991041250257e-06, 1.280761588017388e-06,-4.466710328588646e-07, 4.336940386156499e-08};

      elseif order==16 then
        F := {2.255119743015102e-03, 2.468348151398869e-02, 1.167180741866501e-01, 3.042770443552129e-01, 4.506789844485660e-01,
              3.113322263346599e-01,-6.346360393538239e-02,-2.312686847515921e-01,-1.974115428944359e-02, 1.493343718134866e-01,
              1.933248589902054e-02,-9.361266861388025e-02,-4.412150270924176e-03, 5.368654216335639e-02,-5.366215238458916e-03,
              -2.608403615494328e-02, 7.281544962441208e-03, 9.895088855116517e-03,-4.942686698396113e-03,-2.576894832907076e-03,
              2.211846544558897e-03, 2.884267211839090e-04,-6.654028602161817e-04, 8.078095351238076e-05, 1.235932463637120e-04,
              -4.315894560627099e-05,-9.861077109745385e-06, 8.016192860045321e-06,-7.379163728004056e-07,-5.206891647322203e-07,
              1.632556884112332e-07,-1.491528356269543e-08};

      elseif order==17 then
        F := {1.585196932541634e-03, 1.837444809958346e-02, 9.278294792150479e-02, 2.618775084651301e-01, 4.320398502317169e-01,
              3.665045915597694e-01, 1.931460080002458e-02,-2.321578275719257e-01,-8.951954328819546e-02, 1.395196558810997e-01,
              7.151359021601624e-02,-8.967223561681700e-02,-4.036972996866180e-02, 5.735059315825757e-02, 1.577720421595056e-02,
              -3.317917437455878e-02,-2.312914840193763e-03, 1.607513687385160e-02,-2.151718850854002e-03,-6.083184145000919e-03,
              2.098690587125012e-03, 1.627197831631058e-03,-1.016003058537950e-03,-2.320247296013995e-04, 3.107489840717279e-04,
              -1.810908214139460e-05,-5.801671982744220e-05, 1.639555327116753e-05, 4.943101361106448e-06,-3.186182481274301e-06,
              2.133022685008432e-07, 2.091410386665782e-07,-5.956631070522008e-08, 5.138893560284761e-09};

      elseif order==18 then
        F := {1.114619644715674e-03, 1.363905158130713e-02, 7.324810663591516e-02, 2.225116133165308e-01, 4.043426134369540e-01,
              4.043248276663823e-01, 1.041024608215512e-01,-2.076447635281544e-01,-1.530751364330310e-01, 1.057364881402392e-01,
              1.181443292646286e-01,-6.528850140308115e-02,-7.548523752044420e-02, 4.588219059556805e-02, 4.034132415080043e-02,
              -3.148473688013113e-02,-1.678191401061436e-02, 1.885903701944815e-02, 4.428021425306536e-03,-9.228790681933342e-03,
              8.388410140353035e-05, 3.495471785190439e-03,-7.910634551505644e-04,-9.479447333809070e-04, 4.443923461716852e-04,
              1.510249707643067e-04,-1.404657384127381e-04,-1.086057413370059e-07, 2.645454675509065e-05,-6.024975834038189e-06,
              -2.356528439241738e-06, 1.250668944698616e-06,-5.438805633427158e-08,-8.316274135127112e-08, 2.169994649112698e-08,
              -1.773377459869623e-09};

      elseif order==19 then
        F := {7.839479076417716e-04, 1.009826155732177e-02, 5.747230505201813e-02, 1.869508529511459e-01, 3.708325188059300e-01,
              4.254693669586074e-01, 1.844805901967784e-01,-1.612849715799578e-01,-2.021184348395035e-01, 5.278712614153072e-02,
              1.501539434750123e-01,-2.370118827452947e-02,-1.009647332181425e-01, 1.950508138216296e-02, 6.145235618457190e-02,
              -1.873920386214199e-02,-3.229655512599267e-02, 1.529031257010945e-02, 1.370058271581248e-02,-9.891284492494852e-03,
              -4.148540529617267e-03, 4.978560207890395e-03, 5.437328418525390e-04,-1.900386103064720e-03, 2.416952167285567e-04,
              5.202909518628883e-04,-1.843258632319724e-04,-8.810606477646472e-05, 6.159798420115889e-05, 3.610452213811971e-06,
              -1.176638149984875e-05, 2.129073285963080e-06, 1.083239135480841e-06,-4.852701063231752e-07, 1.023245949054923e-08,
              3.278810145193826e-08,-7.894154721312056e-09, 6.128387585570079e-10};

      elseif order==20 then
        F := {5.515104892333055e-04, 7.459548476695399e-03, 4.484738525098346e-02, 1.555225599602786e-01, 3.342466780734575e-01,
              4.316839091208648e-01, 2.556207268524259e-01,-9.843781145578548e-02,-2.310731625885285e-01,-1.182783757281863e-02,
              1.614261501185060e-01, 2.817837950211079e-02,-1.099259368196373e-01,-1.747743622039748e-02, 7.233116828764272e-02,
              3.982599946733219e-03,-4.364468087822997e-02, 4.154027347113409e-03, 2.283551819224222e-02,-6.214991255579771e-03,
              -9.765516683348944e-03, 4.752908245917196e-03, 3.125795498301260e-03,-2.532498877799901e-03,-5.880032514060865e-04,
              9.846883500163356e-04,-3.782851463836977e-05,-2.723101792733945e-04, 7.179459484541041e-05, 4.790139911371629e-05,
              -2.623780652455027e-05,-3.094401000403553e-06, 5.120335768453775e-06,-7.155878270028505e-07,-4.841616414350393e-07,
              1.862465681522111e-07, 1.424340762347018e-10,-1.283287967659915e-08, 2.868114946327296e-09,-2.120497617473854e-10};

      else
        F := {0};  // error

      end if;

      //----------------------------------------------------
      // filter bank
      (lod, hid, lor, hir) := General.filterBank(F);
    annotation (Documentation(info="<html>
<p>
This function generate the Daubechies wavelet filter banks for the orders up to 10. 
If a wrong order number is given, an empty vector is returned.
</p>
</html>"));
    end wavDaubechies;

    function wavSymlets "Symlets wavelet filters"
      input Integer order=2 "Order of the wavelet. order<=8";
      output Real F[:] "Scaling filter";
      output Real lod[:] "High pass filter for decomposition";
      output Real hid[:] "Low pass filter for decomposition";
      output Real lor[:] "High pass filter for reconstruction";
      output Real hir[:] "Low pass filter for reconstruction";

    algorithm
      if order==1 then
        F := {5.000000000000000e-01, 5.000000000000000e-01};

      elseif order==2 then
        F := {3.415063509462200e-01, 5.915063509458700e-01, 1.584936490537800e-01,-9.150635094586999e-02};

      elseif order==3 then
        F := {2.352336038927000e-01, 5.705584579173100e-01, 3.251825002637100e-01,-9.546720778426000e-02,-6.041610415535000e-02,
              2.490874986589000e-02};

      elseif order==4 then
        F := {2.278517294800000e-02,-8.912350720850001e-03,-7.015881208950001e-02, 2.106172671020000e-01, 5.683291217050001e-01,
              3.518695343280000e-01,-2.095548256255000e-02,-5.357445070900000e-02};

      elseif order==5 then
        F := {1.381607647893000e-02,-1.492124993438000e-02,-1.239756813067500e-01, 1.173946156807000e-02, 4.482908241909200e-01,
              5.115264834460500e-01, 1.409953484272900e-01,-2.767209305836000e-02, 2.087343221079000e-02, 1.932739797744000e-02};

      elseif order==6 then
        F := {-5.515933754690000e-03, 1.249961046390000e-03, 3.162528132994000e-02,-1.489187564922000e-02,-5.136248493090000e-02,
              2.389521856660500e-01, 5.569463919639600e-01, 3.472289864783500e-01,-3.416156079324000e-02,-8.343160770584000e-02,
              2.468306185920000e-03, 1.089235016328000e-02};

      elseif order==7 then
        F := {7.260697381010000e-03, 2.835671342880000e-03,-7.623193594814000e-02,-9.902835340368001e-02, 2.040919698628700e-01,
              5.428913549059901e-01, 3.790813009826900e-01, 1.233282974432000e-02,-3.503914561106000e-02, 4.800738396784000e-02,
              2.157772629104000e-02,-8.935215825569999e-03,-7.406129573000000e-04, 1.896329267100000e-03};

      elseif order==8 then
        F := {1.336396696400000e-03,-2.141971501200000e-04,-1.057284326418000e-02, 2.693194376880000e-03, 3.474523295559000e-02,
              -1.924676063167000e-02,-3.673125438038000e-02, 2.576993351865400e-01, 5.495533152690100e-01, 3.403726735943900e-01,
              -4.332680770282000e-02,-1.013243276428200e-01, 5.379305875240000e-03, 2.241181152181000e-02,-3.833454481100000e-04,
              -2.391729255750000e-03};

      elseif order==9 then
        F := {7.562436546810931e-04,-3.345707545656016e-04,-7.257789276472252e-03, 6.264448120928858e-03, 4.389562577713971e-02,
              -1.289322296471236e-02,-1.354468917522294e-01, 2.494141547906326e-02, 4.365242036747413e-01, 5.076298954167544e-01,
              1.688294618011267e-01,-3.858608054872847e-02, 4.125704643553122e-04, 2.137221680122829e-02,-8.151675612794035e-03,
              -9.384698418212253e-03, 4.382512694515118e-04, 9.905968682437781e-04};

      elseif order==10 then
        F := {-3.247949483908787e-04, 4.033060149896469e-05, 3.247864189340875e-03,-5.687676553368694e-04,-1.439311597192915e-02,
              4.076408391889564e-03, 3.535178378114447e-02,-2.262038615210758e-02,-2.512827017030121e-02, 2.714065055514034e-01,
              5.441257653687310e-01, 3.335356692145733e-01,-5.012010750646451e-02,-1.127794861599803e-01, 8.209434708171549e-03,
              3.247546230148205e-02,-1.036181960273342e-03,-6.110321317044863e-03, 6.762250997113808e-05, 5.445852236221847e-04};

      elseif order==11 then
        F := {3.459616166040179e-04, 7.816011710229956e-05,-4.518132081644625e-03,-1.416668566484035e-03, 3.040572642324027e-02,
              2.493736483851123e-02,-1.022492978094926e-01,-1.447127929903163e-01, 1.680721465083360e-01, 5.164308761562493e-01,
              4.044813267454531e-01, 6.872964384234002e-02,-1.614512237052043e-02, 4.948106953048136e-02, 2.618940799627895e-02,
              -1.702772638911342e-02,-6.970612565931848e-03, 4.605029854079007e-03, 4.160282089289136e-04,-1.226382148667494e-03,
              -2.743267125161877e-05, 1.214257558180884e-04};

      elseif order==12 then
        F := {-1.266191979340108e-04,-1.283970069706689e-05, 1.661911380795638e-03, 2.175398429881296e-04,-1.031657228965601e-02,
              -1.841582559115546e-03, 4.087372726684106e-02, 1.081996455811028e-02,-1.204702753330783e-01,-5.538952842799313e-02,
              3.272073211084249e-01, 5.398612473370130e-01, 2.820549759972962e-01,-1.567111697979455e-02,-2.534895131171365e-02,
              3.477502946382079e-02, 5.341329494112453e-03,-1.712663724874185e-02,-9.962492807415999e-04, 5.243172399797933e-03,
              1.274306051645929e-04,-9.544214477102665e-04,-8.028439511261475e-06, 7.917276232218168e-05};

      elseif order==13 then
        F := {4.980143648796826e-05, 2.609603980976824e-05,-5.100816484363065e-04, 2.922197961024179e-04, 4.012727576024678e-03,
              -1.055319588192468e-03,-1.467224390826674e-02, 1.245801719726598e-02, 6.570862659752902e-02, 6.236510457172705e-03,
              -9.934149753018953e-02, 7.794453819010856e-02, 4.557758467678603e-01, 4.919618712990166e-01, 1.397984180263695e-01,
              -8.793753932244956e-02,-4.225007403951227e-02, 9.802265941068396e-03,-1.217046928712762e-02,-1.429541384079710e-02,
              3.745091886855507e-03, 5.321845009984820e-03,-1.208748544621021e-04,-8.033181615241826e-04,-2.527102313226435e-05,
              4.822698243418468e-05};

      elseif order==14 then
        F := {3.155038190739735e-05, 1.366767897001262e-05,-4.283371327948703e-04,-5.177026689230589e-05, 3.205086977345326e-03,
              7.163995694425472e-04,-1.374567093743074e-02,-1.672342070294283e-03, 4.937558104354086e-02, 1.831306686634989e-02,
              -1.131352543940332e-01,-4.109126433503728e-02, 3.361131410992960e-01, 5.373843542315079e-01, 2.780354625522295e-01,
              -2.497367657535566e-02,-4.075374461460186e-02, 2.646919062210359e-02, 3.026785071847736e-03,-2.064484356596472e-02,
              -1.947212828734811e-03, 7.097721295243721e-03, 2.591380703810256e-04,-1.823940736083744e-03,-4.445256818865413e-05,
              2.817365662280872e-04, 7.927279236449717e-06,-1.829928021760973e-05};

      elseif order==15 then
        F := {2.026618135187328e-05, 1.535686739869793e-05,-2.843761001190123e-04,-7.647671084706012e-05, 2.461459025484652e-03,
              1.079142725531503e-03,-1.214190928336907e-02,-6.183499521675345e-03, 4.806192703374441e-02, 4.836137333158847e-02,
              -9.479211769919799e-02,-1.390358316436560e-01, 1.725076833677631e-01, 5.104201012081355e-01, 4.091605614672253e-01,
              7.886623216607026e-02,-2.904983216720763e-02, 2.880433392844821e-02, 1.551225593035679e-02,-2.748999013390402e-02,
              -1.372141517182688e-02, 7.127620153064661e-03, 2.420745230733134e-03,-2.538630333419648e-03,-1.890212720229747e-04,
              7.570053405981886e-04, 3.897752738499745e-05,-1.136050951918606e-04,-5.204070300795476e-06, 6.867717858445659e-06};

      elseif order==16 then
        F := {-7.635326369100278e-06,-3.815889850668213e-06, 1.169956222972973e-04, 2.585601341018369e-05,-9.466184191064203e-04,
              -1.570600665427034e-04, 4.905737941845218e-03, 9.615554387752017e-04,-1.764426442383797e-02,-2.482139304669221e-03,
              5.518109497727973e-02, 2.172312575746245e-02,-1.128487213186205e-01,-3.821247570060701e-02, 3.361181215191123e-01,
              5.349439490642222e-01, 2.808083193279862e-01,-2.444767136802021e-02,-4.736416822227104e-02, 2.286294833450929e-02,
              3.443096951098992e-03,-2.195651609472216e-02,-2.210781494031773e-03, 8.956731852153130e-03, 5.078525605627896e-04,
              -2.744219371007138e-03,-7.668263492545918e-05, 6.027057959864668e-04, 1.985455582931609e-05,-7.737974108497493e-05,
              -2.201616849456765e-06, 4.405279985266950e-06};

      elseif order==17 then
        F := {2.680820842903412e-06,-1.734332358172868e-06,-5.379049252527110e-05, 1.782470046364011e-05, 5.089945983874696e-04,
              4.129533915321227e-05,-2.780573871173983e-03,-1.347326698416031e-03, 8.765994540344152e-03, 7.037821435399957e-03,
              -1.275542656487228e-02,-5.134751174859099e-03, 1.142600322601135e-02,-6.086129917757593e-02,-1.096552949818622e-01,
              1.276607645300947e-01, 4.818854899122590e-01, 4.601261411410810e-01, 1.006908392078749e-01,-8.383948207398410e-02,
              1.221256723174161e-02, 7.407269851545650e-02, 1.266000602076373e-02,-2.354056302252066e-02,-3.407698053101591e-03,
              7.412152741227024e-03, 6.058279265480015e-04,-1.938657674219570e-03,-9.803491238421889e-05, 3.365802683046743e-04,
              -9.550455291405712e-06,-4.450319785498648e-05, 1.965846437781149e-06, 3.038680607852918e-06};

      elseif order==18 then
        F := {-1.069960796228225e-06, 5.548877669031676e-07, 2.090026457626608e-05,-6.971235668653580e-06,-1.879702735956231e-04,
              3.352827779049330e-05, 1.009809525994678e-03,-1.334849589823968e-04,-3.705090816839179e-03, 7.691800010747126e-04,
              1.061533897250874e-02,-2.305694335338614e-03,-2.242425442343012e-02, 4.439177165357097e-03, 2.017347152970506e-02,
              -5.218391992277767e-02,-2.296723363564539e-02, 2.838919603368541e-01, 5.328962754611726e-01, 3.351467363307968e-01,
              -3.679017112570010e-02,-1.130933494899401e-01, 2.403856673886823e-02, 5.955248359256860e-02,-3.590041344561848e-03,
              -2.144307754874283e-02, 1.161766822651304e-03, 6.719044876650879e-03,-2.909893668402564e-04,-1.636154450781406e-03,
              4.964790073697408e-05, 2.801333666356325e-04,-9.914338927982427e-06,-3.199428931879830e-05, 9.580701229228587e-07,
              1.847396055208554e-06};

      elseif order==19 then
        F := {1.238099284413472e-06, 1.458278380500544e-06,-1.990586104598978e-05,-1.189451683876792e-05, 1.953161702093965e-04,
              9.143433492874425e-05,-1.205588962334967e-03,-4.369370682946555e-04, 5.842283779057898e-03, 3.054243001324849e-03,
              -1.959385602616820e-02,-1.195592752929102e-02, 5.944835950907114e-02, 6.620700430926962e-02,-8.219531561490194e-02,
              -1.248728416684554e-01, 1.826217596092635e-01, 5.088025916704174e-01, 4.088102113629857e-01, 7.709290024687236e-02,
              -4.774742643691535e-02, 6.331852141968008e-03, 4.960759849055843e-03,-3.297662020450615e-02,-1.601737812555945e-02,
              1.117047645102727e-02, 5.634536772047775e-03,-3.621945891862454e-03,-8.207411441384112e-04, 1.500074023194426e-03,
              1.125417348133064e-04,-4.495533997908204e-04,-3.261204476796456e-05, 8.169857538311178e-05, 6.274379210031963e-06,
              -8.400795032381557e-06,-4.570491667968555e-07, 3.880413053636296e-07};

      elseif order==20 then
        F := {-4.475370066899557e-07,-2.302836523067074e-07, 8.688399761964284e-06, 3.199956732525194e-06,-8.300820915085660e-05,
              -1.882003613225947e-05, 5.286407086683141e-04, 8.870012323726663e-05,-2.454825703404990e-03,-4.321316116626793e-04,
              8.596326094449493e-03, 1.370795132503538e-03,-2.501272629501434e-02,-4.839228068872543e-03, 6.287570024610344e-02,
              2.563329374068724e-02,-1.135460037275814e-01,-3.612491372832120e-02, 3.337483727097380e-01, 5.311522590562131e-01,
              2.869661663249030e-01,-2.108547794998238e-02,-5.585743697906758e-02, 1.808733149590397e-02, 5.743989855459407e-03,
              -2.236538949158362e-02,-2.343251028255893e-03, 1.202367837243984e-02, 1.006274722107307e-03,-4.671561619256600e-03,
              -2.158534184425312e-04, 1.477142324081319e-03, 5.102476302454971e-05,-3.498277097257839e-04,-1.363593414806457e-05,
              5.651881758718569e-05, 2.139468990693039e-06,-5.599834157312889e-06,-1.344611337171996e-07, 2.613139608698177e-07};

      end if;

      // filter bank
      (lod, hid, lor, hir) := General.filterBank(F);
    annotation (Documentation(info="<html>
<p>
This function generate the Symlets wavelet filter banks for the orders up to 8.
If a wrong order number is given, an empty vector is returned.
</p>
</html>"));
    end wavSymlets;

    function wavCoiflets "Coiflets wavelet filters"
      input Integer order=2 "Order of the wavelet. order<=5";
      output Real F[:] "Scaling filter";
      output Real lod[:] "High pass filter for decomposition";
      output Real hid[:] "Low pass filter for decomposition";
      output Real lor[:] "High pass filter for reconstruction";
      output Real hir[:] "Low pass filter for reconstruction";

    algorithm
      if order==1 then
        F :={-0.051429728471000, 0.238929728471000, 0.602859456942000, 0.272140543058000,-0.051429728471000,
             -0.011070271529000};

      elseif order==2 then
        F :={0.011587596739000,-0.029320137980000,-0.047639590310000, 0.273021046535000, 0.574682393857000,
             0.294867193696000,-0.054085607092000,-0.042026480461000, 0.016744410163000, 0.003967883613000,
             -0.001289203356000,-0.000509505399000};

      elseif order==3 then
        F :={-0.002682418671000, 0.005503126709000, 0.016583560479000,-0.046507764479000,-0.043220763560000,
             0.286503335274000, 0.561285256870000, 0.302983571773000,-0.050770140755000,-0.058196250762000,
             0.024434094321000, 0.011229240962000,-0.006369601011000,-0.001820458916000, 0.000790205101000,
             0.000329665174000,-0.000050192775000,-0.000024465734000};

      elseif order==4 then
        F :={0.000630961046000,-0.001152224852000,-0.005194524026000, 0.011362459244000, 0.018867235378000,
             -0.057464234429000,-0.039652648517000, 0.293667390895000, 0.553126452562000, 0.307157326198000,
             -0.047112738865000,-0.068038127051000, 0.027813640153000, 0.017735837438000,-0.010756318517000,
             -0.004001012886000, 0.002652665946000, 0.000895594529000,-0.000416500571000,-0.000183829769000,
             0.000044080354000, 0.000022082857000,-0.000002304942000,-0.000001262175000};

      elseif order==5 then
        F :={-0.000149963800000, 0.000253561200000, 0.001540245700000,-0.002941110800000,-0.007163781900000,
             0.016552066400000, 0.019917804300000,-0.064997262800000,-0.036800073600000, 0.298092323500000,
             0.547505429400000, 0.309706849000000,-0.043866050800000,-0.074652238900000, 0.029195879500000,
             0.023110777000000,-0.013973687900000,-0.006480090000000, 0.004783001400000, 0.001720654700000,
             -0.001175822200000,-0.000451227000000, 0.000213729800000, 0.000099377600000,-0.000029232100000,
             -0.000015072000000, 0.000002640800000, 0.000001459300000,-0.000000118400000,-0.000000067300000};

      end if;

      // filter bank
      (lod, hid, lor, hir) := General.filterBank(F);
    annotation (Documentation(info="<html>
<p>
This function generate the Coiflets wavelet filter banks for the orders up to 5.
If a wrong order number is given, an empty vector is returned.
</p>
</html>"));
    end wavCoiflets;

    function wavBiorSpline "Biorthogonal spline wavelet filters"
      import Modelica.Math.Vectors.reverse;
      input Integer Nd=2 "Wavelet order for decomposition";
      input Integer Nr=2 "Wavelet order for reconstruction";
      output Real Fd[:] "Scaling filter for decomposition";
      output Real lodFd[:] "High pass filter for decomposition using filter Fd";
      output Real hidFd[:] "Low pass filter for decomposition using filter Fd";
      output Real lorFd[:]
        "High pass filter for reconstruction using filter Fd";
      output Real hirFd[:] "Low pass filter for reconstruction using filter Fd";
      output Real Fr[:] "Scaling filter for reconstruction";
      output Real lodFr[:] "High pass filter for decomposition using filter Fr";
      output Real hidFr[:] "Low pass filter for decomposition using filter Fr";
      output Real lorFr[:]
        "High pass filter for reconstruction using filter Fr";
      output Real hirFr[:] "Low pass filter for reconstruction using filter Fr";

    protected
      Integer nFd "Fd length";
      Integer nFr "Fr length";
      Integer nF "filter length";

    algorithm
      //------------------------------------
      if Nr==1 then
        Fr :={1/2};

        if Nd==1 then
          Fd :={1/2};
        elseif Nd==3 then
          Fd :={-1/16, 1/16, 1/2};
        elseif Nd==5 then
          Fd :={3/256, -3/256, -11/128, 11/128, 1/2};
        end if;

        Fr :=cat(1, Fr, reverse(Fr));
        Fd :=cat(1, Fd, reverse(Fd));

      //------------------------------------
      elseif Nr==2 then
        Fr :={1/4,1/2,1/4};

        if Nd == 2 then
           Fd :={-1/8,1/4};
           Fd :=cat(
            1,
            Fd,
            {3/4},
            reverse(Fd));
        elseif Nd == 4 then
           Fd :={3/128,-3/64,-1/8,19/64};
           Fd :=cat(
            1,
            Fd,
            {45/64},
            reverse(Fd));
        elseif Nd == 6 then
           Fd :={-5/1024,5/512,17/512,-39/512,-123/1024,81/256};
           Fd :=cat(
            1,
            Fd,
            {175/256},
            reverse(Fd));
        elseif Nd == 8 then
           Fd :={35,-70,-300,670,1228,-3126,-3796,10718};
           Fd :=cat(
            1,
            Fd,
            {22050},
            reverse(Fd))/(2^15);
        end if;

      //------------------------------------
      elseif Nr==3 then
        Fr :={1,3}/8;

        if Nd == 1 then
           Fd :={-1,3}/4;
        elseif Nd == 3 then
           Fd :={3,-9,-7,45}/64;
        elseif Nd == 5 then
           Fd :={-5,15,19,-97,-26,350}/512;
        elseif Nd == 7 then
           Fd :={35,-105,-195,865,363,-3489,-307,11025}/(2^14);
        elseif Nd == 9 then
           Fd :={-63,189,469,-1911,-1308,9188,1140,-29676,190,87318}/(2^17);
        end if;

        Fr :=cat(
          1,
          Fr,
          reverse(Fr));
        Fd :=cat(
          1,
          Fd,
          reverse(Fd));

      //------------------------------------
      elseif Nr==4 and Nd==4 then
           Fr :={-0.045635881557,-0.028771763114,0.295635881557};
           Fr :=cat(
          1,
          Fr,
          {0.557543526229},
          reverse(Fr));

           Fd :={0.026748757411,-0.016864118443,-0.078223266529,0.266864118443};
           Fd :=cat(
          1,
          Fd,
          {0.602949018236},
          reverse(Fd));

      //------------------------------------
      elseif Nr==5 and Nd == 5 then
           Fr :={0.009515330511,-0.001905629356,-0.096666153049,-0.066117805605,0.337150822538};
           Fr :=cat(
          1,
          Fr,
          {0.636046869922},
          reverse(Fr));

           Fd :={0.028063009296,0.005620161515,-0.038511714155,0.244379838485};
           Fd :=cat(
          1,
          Fd,
          {0.520897409718},
          reverse(Fd));

      //------------------------------------
      elseif Nr==6 and Nd == 8 then
    /*       Fr :={-0.01020092218704,-0.01023007081937,0.05566486077996,0.02854447171515,
      -0.29546393859292};
       Fr :=cat(
      1,
      Fr,
      {-0.53662880179157},
      reverse(Fr));

       Fd :={0.00134974786501,-0.00135360470301,-0.01201419666708,0.00843901203981,
      0.03516647330654,-0.05463331368252,-0.06650990062484,0.29754790634571};
       Fd :=cat(
      1,
      Fd,
      {0.58401575224075},
      reverse(Fd));
*/

          Fr:={-0.010200922187040,-0.010230070819370,0.055664860779960,0.028544471715150,-0.295463938592920,
          -0.536628801791570,-0.295463938592920,0.028544471715150,0.055664860779960,-0.010230070819370,-0.010200922187040};

          Fd:={0.001349747865010,-0.001353604703010,-0.012014196667080,0.008439012039810,0.035166473306540,
          -0.054633313682520,-0.066509900624840,0.297547906345710,0.584015752240750,0.297547906345710,-0.066509900624840,
          -0.054633313682520,0.035166473306540,0.008439012039810,-0.01201419666708,-0.001353604703010,0.001349747865010};
      end if;

      // make both filters the same length
      nFr :=size(Fr, 1);
      nFd :=size(Fd, 1);
      nF :=max(nFr, nFd);
      if mod(nF,2)>0 then
         nF :=nF + 1;
      end if;
      Fr :=cat(
        1,
        fill(0., integer(floor((nF - nFr)/2))),
        Fr,
        fill(0., integer(ceil((nF - nFr)/2))));
      Fd :=cat(1,
        fill(0., integer(floor((nF - nFd)/2))),
        Fd,
        fill(0., integer(ceil((nF - nFd)/2))));

      // filter banks
      (lodFd, hidFd, lorFd, hirFd) := General.filterBank(Fd);
      (lodFr, hidFr, lorFr, hirFr) := General.filterBank(Fr);

    annotation (Documentation(info="<html>
<p>
This function generates the biorthogonal spline wavelet filters, where Nr is the order of the reconstruction filters, 
and Nd the order of the decomposition filters.
</p>
<p>
Usually, lodFd and hidFr are used for decomposition, and lorFr and hirFd are used for reconstruction.
</p>
<h4><font color=\"#008000\">Available orders of biorthogonal wavelet filters</font></h4>

<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><i>Nd</i></td>
      <td><i>Nr</i></td>
      </tr>
  
  <tr><td>1, 3, 5</td> <td>1</td>
      </tr>
  
  <tr><td>2, 4, 6, 8 </td> <td>2</td>
      </tr>
  
  <tr><td>1, 3, 5, 7, 9</td> <td>3</td>
      </tr>
  
  <tr><td>4 </td> <td>4</td>
      </tr>
  
  <tr><td>5 </td> <td>5</td>
      </tr>
  
  <tr><td>8 </td> <td>6</td>
      </tr>
  
</table>

</html>"));
    end wavBiorSpline;

    function wavRevBiorSpline "Reverse biorthogonal spline wavelet filters"
      input Integer Nd=2 "Wavelet order for decomposition";
      input Integer Nr=2 "Wavelet order for reconstruction";
      output Real Fd[:] "Scaling filter for decomposition";
      output Real lodFd[:] "High pass filter for decomposition using filter Fd";
      output Real hidFd[:] "Low pass filter for decomposition using filter Fd";
      output Real lorFd[:]
        "High pass filter for reconstruction using filter Fd";
      output Real hirFd[:] "Low pass filter for reconstruction using filter Fd";
      output Real Fr[:] "Scaling filter for reconstruction";
      output Real lodFr[:] "High pass filter for decomposition using filter Fr";
      output Real hidFr[:] "Low pass filter for decomposition using filter Fr";
      output Real lorFr[:]
        "High pass filter for reconstruction using filter Fr";
      output Real hirFr[:] "Low pass filter for reconstruction using filter Fr";

    algorithm
      (Fr, lodFr, hidFr, lorFr, hirFr, Fd, lodFd, hidFd, lorFd, hirFd):= wavBiorSpline(Nr, Nd);

    annotation (Documentation(info="<html>
<p>
This function generates the biorthogonal spline wavelet filters, where Nr is the order of the reconstruction filters, 
and Nd the order of the decomposition filters.
</p>
<p>
Usually, lodFd and hidFr are used for decomposition, and lorFr and hirFd are used for reconstruction.
</p>

<h4><font color=\"#008000\">Available orders of biorthogonal wavelet filters</font></h4>

<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><i>Nd</i></td>
      <td><i>Nr</i></td>
      </tr>
  
  <tr><td>1</td> <td>1, 3, 5</td>
      </tr>
  
  <tr><td>2 </td> <td>2, 4, 6, 8</td>
      </tr>
  
  <tr><td>3 </td> <td>1, 3, 5, 7, 9</td>
      </tr>
  
  <tr><td>4 </td> <td>4</td>
      </tr>
  
  <tr><td>5 </td> <td>5</td>
      </tr>
  
  <tr><td>6 </td> <td>8</td>
      </tr>
  
</table>

</html>"));
    end wavRevBiorSpline;

    function wavDMeyer "Discrete Meyer wavelet filters"

      input Integer points= 22
        "Filter length (<= 102). The complete length is used if this variable is <= 0";
        output Real F[:] "Scaling filter";
        output Real lod[:] "High pass filter for decomposition";
        output Real hid[:] "Low pass filter for decomposition";
        output Real lor[:] "High pass filter for reconstruction";
        output Real hir[:] "Low pass filter for reconstruction";

    protected
      Integer id1 "start index";
      Integer id2 "end index";
      Real F0[:] = {-1.06754713000000e-06, 9.04223910000000e-07, 3.17904740000000e-07, -1.48249686000000e-06, 1.21850207000000e-06,
             4.93618310000000e-07, -2.03604729000000e-06, 1.68513902000000e-06, 6.94742880000000e-07, -2.98242491000000e-06,
             2.37128175000000e-06, 1.18420622000000e-06, -4.26703335000000e-06, 3.42066573000000e-06, 1.69867277000000e-06,
             -6.75732600000000e-06, 5.10285152000000e-06, 3.42881336000000e-06, -1.00458073700000e-05, 7.42738297000000e-06,
             4.37527643000000e-06, -1.72802656000000e-05, 1.42173515200000e-05, 1.06020135900000e-05, -3.28300673700000e-05,
             2.28687423700000e-05, 2.64526068300000e-05, -7.26756723600000e-05, 1.72972015000000e-05, 0.000105863355880000,
             -5.34521877000000e-05, -9.89334554300000e-05, -6.61235476200000e-05, 0.000113978321900000, 0.000607757935360000,
             -0.000408838764160000, -0.00191072028190000, 0.00155193926157000, 0.00427481806225000, -0.00451609544333000,
             -0.00780973483285000, 0.0107840153442700, 0.0123063973649900, -0.0226939113793400, -0.0171980843830200,
             0.0450195435858100, 0.0216524716331900, -0.0938306002586500, -0.0247828615296900, 0.314022352386430,
             0.525910951416260, 0.314022352386430, -0.0247828615296900, -0.0938306002586500, 0.0216524716331900,
             0.0450195435858100, -0.0171980843830200, -0.0226939113793400, 0.0123063973649900, 0.0107840153442700,
             -0.00780973483285000, -0.00451609544333000, 0.00427481806225000, 0.00155193926157000, -0.00191072028190000,
             -0.000408838764160000, 0.000607757935360000, 0.000113978321900000, -6.61235476200000e-05, -9.89334554300000e-05,
             -5.34521877000000e-05, 0.000105863355880000, 1.72972015000000e-05, -7.26756723600000e-05, 2.64526068300000e-05,
             2.28687423700000e-05, -3.28300673700000e-05, 1.06020135900000e-05, 1.42173515200000e-05, -1.72802656000000e-05,
             4.37527643000000e-06, 7.42738297000000e-06, -1.00458073700000e-05, 3.42881336000000e-06, 5.10285152000000e-06,
             -6.75732600000000e-06, 1.69867277000000e-06, 3.42066573000000e-06, -4.26703335000000e-06, 1.18420622000000e-06,
             2.37128175000000e-06, -2.98242491000000e-06, 6.94742880000000e-07, 1.68513902000000e-06, -2.03604729000000e-06,
             4.93618310000000e-07, 1.21850207000000e-06, -1.48249686000000e-06, 3.17904740000000e-07, 9.04223910000000e-07,
             -1.06754713000000e-06, 0};

    algorithm
        //----------------------------------------------------
        // get index
        if points <=0 then
          id1 := 1;
          id2 := 102;
        else
          id1 := integer(51 + floor(-points/2.0)+1);
          id2 := integer(51 + floor(points/2.0));
          if id1<1 then
            id1:=1;
          end if;
          if id2 > 102 then
            id2 := 102;
          end if;
        end if;

        // filter bank
        F := F0[id1:id2];
        (lod, hid, lor, hir) := General.filterBank(F);

    annotation (Documentation(info="<html>
<p>
This function returns the filter bank of discrete Meyer wavelet, which is an approximation of continuous Meyer wavelet. 
The filters with the length of 102 points usually provide sufficient accuracy. However, such long filters might cause 
long calculation time. Therefore, this function provides the possibility to return shorter filters. 
The trade-off is lower accuracy. Please note that some lengths will cause significantly large errors. The default length is 22. 
</p>

</html>"));
    end wavDMeyer;

    function wavMeyer
      "Meyer wavelet function. The fft_c.c file must be accessible."
      import Modelica.ComplexMath.*;

        input Real low=-8 "Lower boundary of the function";
        input Real high=8 "Higher boundary of the function";
        input Integer points=32
        "Number of data points of the returned functions";

        output Real x[points] "Regular grid";
        output Real phi[points] "Scaling function";
        output Real psi[points] "Wavelet function";
    protected
        Real pi = Modelica.Constants.pi;
        Real lint;
        Real[:] x1;
        Real[:] xa;
        Integer n_xa;
        Complex[:] phihat;
        Complex[:] psihat;
        Real tmp;
    algorithm
        lint := (high-low)/2/pi;
        x1 := {(i/(2*lint)) for i in (-points):2:(points-2)};
        xa := abs(x1);

        n_xa := size(xa, 1);

        // scaling function phi
        phihat := {Complex(0) for i in 1:points};
        for i in 1:n_xa loop
            if (xa[i] < 2*pi/3) then
                phihat[i] := Complex(1);
            end if;
            if (xa[i] >= 2*pi/3 and xa[i] < 4*pi/3) then
                tmp := 3/2/pi*xa[i] - 1;
                tmp := 35*tmp^4 - 84*tmp^5 + 70*tmp^6 - 20*tmp^7;
                phihat[i] := Complex(Modelica.Math.cos(pi/2*tmp));
            end if;
        end for;

        (phi, x) := Wavelet.General.nStdIfft(phihat, low, high);

        // wavelet function psi
        psihat := {Complex(0) for i in 1:points};
        for i in 1:n_xa loop
            if (xa[i] >= 2*pi/3 and xa[i] < 4*pi/3) then
                tmp := 3/2/pi*xa[i] - 1;
                tmp := 35*tmp^4 - 84*tmp^5 + 70*tmp^6 - 20*tmp^7;
                psihat[i] := exp(j*x1[i]/2) * Modelica.Math.sin(pi/2*tmp);
            end if;
            if (xa[i] >= 4*pi/3 and xa[i] < 8*pi/3) then
                tmp := 3/4/pi*xa[i] - 1;
                tmp := 35*tmp^4 - 84*tmp^5 + 70*tmp^6 - 20*tmp^7;
                psihat[i] := exp(j*x1[i]/2) * Modelica.Math.cos(pi/2*tmp);
            end if;
        end for;

        (psi, x) := Wavelet.General.nStdIfft(psihat, low, high);
    annotation (Documentation(info="<html>
<p>
This function generates the continuous Meyer wavelet and scaling functions. Since Meyer wavelet has theoretically 
unlimited support, the function range has to be defined with [low, high].
</p><p>
For generating Meyer wavelets, an external C-function saved in fft_c.c file has to be called. This C-file has to 
be available in the current working directory or the searching path of Modelica must include the directory of fft_c.c file.
</p>
</html>"));
    end wavMeyer;

    function wavGaussian "Gaussian wavelet function"
      import Modelica.Constants.*;
       input Real low=-5 "Lower boundary of the function";
       input Real high=5 "Higher boundary of the function";
       input Integer order=1 "Order of the derivative, <=8";
       input Integer points=32
        "Number of data points of the returned functions";
       output Real x[points] "Regular grid";
       output Real psi[points] "Wavelet function";

    protected
       Integer i;
       Real tem;
       Real k;
       Real X2[points];
       Real F0[points];

    algorithm
       x :=linspace(
       low,
       high,
       points);

       for i in 1:points loop
         X2[i]:=x[i]^2;
         F0[i] :=((2/pi)^(1/4))*Modelica.Math.exp(-X2[i]);
       end for;

       if order==1 then
         for i in 1:points loop
           psi[i]:=-2*x[i]*F0[i];
         end for;

       elseif order==2 then
         for i in 1:points loop
           psi[i]:=2/(3^(1/2))*(-1 + 2*X2[i])*F0[i];
         end for;

        elseif order==3 then
         for i in 1:points loop
           psi[i]:=4/(15^(1/2))*x[i]*(3 - 2*X2[i])*F0[i];
         end for;

         elseif order==4 then
         for i in 1:points loop
           psi[i]:=4/(105^(1/2))*(3 - 12*X2[i] + 4*X2[i]^2)*F0[i];
         end for;

         elseif order==5 then
         for i in 1:points loop
           psi[i]:=8/(3*(105^(1/2)))*x[i]*(-15 + 20*X2[i] - 4*X2[i]^2)*F0[i];
         end for;

         elseif order==6 then
         for i in 1:points loop
           psi[i]:=8/(3*(1155^(1/2)))*(-15 + 90*X2[i] - 60*X2[i]^2 + 8*X2[i]^3)*F0[
            i];
         end for;

         elseif order==7 then
         for i in 1:points loop
           psi[i]:=16/(3*(15015^(1/2)))*x[i]*(105 - 210*X2[i] + 84*X2[i]^2 - 8*X2[i]
            ^3)*F0[i];
         end for;

         elseif order==8 then
         for i in 1:points loop
           psi[i]:=16/(45*(1001^(1/2)))*(105 - 840*X2[i] + 840*X2[i]^2 - 224*X2[i]^3
             + 16*X2[i]^4)*F0[i];
         end for;
         end if;

         i:=rem(order, 4);
         if i==2 then
           psi:=-psi;
         elseif i==3 then
           psi:=-psi;
        end if;
    annotation (Documentation(info="<html>
<p>
This function generates the continuous Gaussian wavelet function with the order up to 8. It has no scaling function. 
Since Gaussian wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavGaussian;

    function wavMexHat "Mexican hat wavelet function"
      import Modelica.Constants.*;
       input Real low=-5 "Lower boundary of the returned function";
       input Real high=5 "Higher boundary of the returned function";
       input Integer points=32
        "Number of the data points of the returned function";
       output Real x[points] "Regular grid";
       output Real psi[points] "The wavelet function";

    protected
       Integer i;
       Real tem;
       Real k;

    algorithm
       x :=linspace(
       low,
       high,
       points);

       for i in 1:points loop
         psi[i]:=x[i]*x[i];
       end for;

       tem:=2/((pi^0.25)*(3^0.5));
       for i in 1:points loop
         psi[i]:=tem*(Modelica.Math.exp(-psi[i]/2))
          *(1 - psi[i]);
       end for;

    annotation (Documentation(info="<html>
<p>
This function generates the continuous Mexican Hat wavelet function. It has no scaling function. 
Since Mexican Hat wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavMexHat;

    function wavMorlet "Morlet wavelet function"
       input Real low=-4 "Lower boundary of the function";
       input Real high=4 "Higher boundary of the function";
       input Integer points=32
        "Integer number, length of the returned functions";
       output Real x[points] "Regular grid";
       output Real psi[points] "Wavelet function";

    protected
       Integer i;
       Real tem;
       Real k;

    algorithm
       x :=linspace(
       low,
       high,
       points);
       for i in 1:points loop
         psi[i]:=Modelica.Math.exp(-(x[i]^2)/2)*Modelica.Math.cos(5*x[i]);
       end for;
    annotation (Documentation(info="<html>
<p>
This function generates the continuous Morlet wavelet function. It has no scaling function. 
Since Morlet wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavMorlet;

    function wavXGaussian "Complex Gaussian wavelet function"
      import Modelica.Constants.*;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.Math.Complex;

       input Real low=-5 "Lower boundary of the function";
       input Real high=5 "Higher boundary of the function";
       input Integer order=2 "Order of the derivative, <=8";
       input Integer points=32
        "Number of data points of the returned functions";
       output Real x[points] "Regular grid";
       output Complex psi[points] "Wavelet function";

    protected
       Integer i;
       Complex j=Complex(0,1);
       Real tem;
       Real k;
       Real X2[points];
       Real F0[points];
       Complex F1[points];
       Complex F2[points];

    algorithm
       x :=linspace(
       low,
       high,
       points);
       // x := low:(high-low)/(points-1):high;

       for i in 1:points loop
         X2[i]:=x[i]^2;
         F0[i] :=Modelica.Math.exp(-X2[i]);
         F1[i]:=Complex(Modelica.Math.cos(-x[i]), Modelica.Math.sin(-x[i]));
         F2[i] :=(F1[i]*F0[i])/(exp(-1/2)*2^(1/2)*pi^(1/2))^(1/2);
       end for;

       if order==1 then
         for i in 1:points loop
           psi[i]:=F2[i]*(-j-2*x[i])*2^(1/2);
         end for;

       elseif order==2 then
         for i in 1:points loop
           psi[i]:=1/3*F2[i]*(-3+4*j*x[i]+4*X2[i])*6^(1/2);
         end for;

        elseif order==3 then
         for i in 1:points loop
           psi[i]:=1/15*F2[i]*(7*j+18*x[i]-12*j*x[i]^2-8*x[i]^3)*30^(1/2);
         end for;

         elseif order==4 then
         for i in 1:points loop
           psi[i]:=1/105*F2[i]*(25-56*j*x[i]-72*x[i]^2+32*j*x[i]^3+16*x[i]^4)*210^(1/2);
         end for;

         elseif order==5 then
         for i in 1:points loop
           psi[i]:=1/315*F2[i]*(-81*j-250*x[i]+280*j*x[i]^2+240*x[i]^3-80*j*x[i]^4-32*x[i]^5)*210^(1/2);
         end for;

         elseif order==6 then
         for i in 1:points loop
           psi[i]:=1/3465*F2[i]*(-331+972*j*x[i]+1500*x[i]^2-1120*j*x[i]^3-720*x[i]^4+192*j*x[i]^5+64*x[i]^6)*2310^(1/2);
         end for;

         elseif order==7 then
         for i in 1:points loop
           psi[i]:=1/45045*F2[i]*(1303*j+4634*x[i]-6804*j*x[i]^2-7000*x[i]^3+3920*j*x[i]^4+2016*x[i]^5-448*j*x[i]^6-128*x[i]^7)*30030^(1/2);
         end for;

         elseif order==8 then
         for i in 1:points loop
           psi[i]:=1/45045*F2[i]*(5937-20848*j*x[i]-37072*x[i]^2+36288*j*x[i]^3+28000*x[i]^4-12544*j*x[i]^5-5376*x[i]^6+1024*j*x[i]^7+256*x[i]^8)*2002^(1/2);
         end for;
         end if;

       tem :=0;
       for i in 1:points loop
         tem:=tem + psi[i].re^2+psi[i].im^2;
       end for;
       for i in 1:points loop
           psi[i]:=psi[i]/((tem*(x[2]-x[1]))^0.5);
       end for;
    annotation (Documentation(info="<html>
<p>
This function generates the complex Gaussian wavelet function with the order up to 8. It has no scaling function. 
Since complex Gaussian wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavXGaussian;

    function wavXMorlet "Complex Morlet wavelet function"
      import Modelica.Constants.pi;
      import Modelica_LinearSystems2.Math.Complex;

       input Real low=-4 "Lower boundary of the function";
       input Real high=4 "Higher boundary of the function";
       input Integer points=32
        "Number of data points of the returned functions";
       input Real fb=1 "Band width";
       input Real fc=1 "Center frequency, e.g. 0.1, 0.5, 1, 1.5, ...";
       output Real x[points] "Regular grid";
       output Complex psi[points] "Wavelet function";

    protected
       Integer i;

    algorithm
      x :=linspace(
      low,
      high,
      points);
    //   x := low:(high-low)/(points-1):high;

      for i in 1:points loop
        psi[i] :=1/sqrt(pi*fb) * Complex(cos(2*pi*fc*x[i]),sin(2*pi*fc*x[i])) * exp(-x[i]*x[i]/fb);
      end for;
    annotation (Documentation(info="<html>
<p>
This function generates the complex Morlet wavelet function. It has no scaling function. 
Since complex Morlet wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavXMorlet;

    function wavXShannon "Complex Shannon wavelet function"
      import Modelica.Constants.*;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.Math.Complex;

       input Real low=-20 "Lower boundary of the function";
       input Real high=20 "Higher boundary of the function";
       input Integer points=32
        "Number of data points of the returned functions";
       input Real fb=1 "Band width, e.g. 1, 2";
       input Real fc=1 "Center frequency, e.g. 0.1, 0.5, 1, 1.5, ...";
       output Real x[points] "Regular grid";
       output Complex psi[points] "Wavelet function";

    protected
       Integer i;

    algorithm
       x :=linspace(
       low,
       high,
       points);
    //    x := low:(high-low)/(points-1):high;

       for i in 1:points loop
         psi[i]:=sqrt(fb) * General.sinc(fb*x[i]) * (Complex(cos(2*pi*fc*x[i]),sin(2*pi*fc*x[i])));
       end for;
    annotation (Documentation(info="<html>
<p>
This function generates the complex Shannon wavelet function. It has no scaling function. 
Since complex Shannon wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavXShannon;

    function wavXFreqBSpline "Complex Frequency B-Spline wavelet"
      import Modelica.Constants.*;
      import Modelica.Utilities.Streams;
      import Modelica_LinearSystems2.Math.Complex;

       input Real low=-20 "Lower boundary of the wavelet function";
       input Real high=20 "Higher boundary of the wavelet function";
       input Integer order=2 "Wavelet order, a natural number, 1, 2, ...";
       input Integer points=32
        "Number of data points of the returned functions";
       input Real fb=1 "Band width";
       input Real fc=1
        "Center frequency, a real number, e.g. 0.1, 0.5, 1, 1.5, ...";
       output Real x[points] "Regular grid";
       output Complex psi[points] "Wavelet function";

    protected
       Integer i;

    algorithm
       x :=linspace(
       low,
       high,
       points);
       // x := low:(high-low)/(points-1):high;

       for i in 1:points loop
         psi[i] := sqrt(fb) * ((General.sinc(fb*x[i]/order))^order) *(Complex(cos(2*pi*fc*x[i]),sin(2*pi*fc*x[i])));
       end for;
    annotation (Documentation(info="<html>
<p>
This function generates the complex Frequency B-Spline wavelet function. It has no scaling function. 
Since this wavelet has theoretically unlimited support, the function range has to be defined with [low, high].
</p>
</html>"));
    end wavXFreqBSpline;

  annotation (Documentation(info="<html>
<p>
This section defines all available wavelets in this library. Since different wavelets have different properties, 
the functions for different wavelets require different input parameters and generate different output data. 
All available wavelets in this library with the parameters are listed in the following table.
</p><p>Parameter, Nd, specifies the wavelet order. Nr is only used for biorthogonal wavelets, where Nd is for 
decomposition and Nr for reconstruction. 
</p><p>For the wavelets with theoretically unlimited support, e.g. Meyer, Gaussian, etc., the function range 
is set as effective support, [-es, es], if 'range' is not specified (==0). 
</p><p>Parameters, fb and fc, are used for some complex wavelets. These are the so-called frequency parameters 
defining the oscillation properties of the wavelets.
</p><p>'phi' stands for scaling function, and 'psi' for wavelet function. A biorthogonal wavelet has two
scaling and wavelet functions. The first function set is for forward transformation (decomposition); 
the second set is for inverse transformation (reconstruction).
</p><p>Symbol, O, means the corresponding parameter is applicable, - not applicable, for a specific wavelet. 
</p><p>
For generating Meyer wavelets, an external C-function saved in fft_c.c file has to be called. This C-file has to 
be available in the current working directory or the searching path of Modelica must include the directory of fft_c.c file.</p>
<h4><font color=\"#008000\">Applicable parameters of the wavelet families</font></h4>

<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Wavelet</b></td>
      <td><b>wavID</b></td>
      <td><b>wavType</b></td>
      <td><b>Nd</b></td>
      <td><b>Nr</b></td>
      <td><b>es</b></td>
      <td><b>range</b></td>
      <td><b>fb</b></td>
      <td><b>fc</b></td>
      <td><b>phi1</b></td>
      <td><b>psi1</b></td>
      <td><b>phi2</b></td>
      <td><b>psi2</b></td>
      </tr>
      
  
  <tr><td>Haar </td>
      <td>1</td> <td>1</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Daubechies </td>
      <td>2</td>  <td>1</td> <td>1..20</td> <td>-</td> <td>-</td> <td>-</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Symlets </td>
      <td>3</td> <td>1</td> <td>1..20</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Coiflets </td>
      <td>4</td> <td>1</td> <td>1..5</td> <td>-</td> <td>-</td> <td>-</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Biorthogonal spline </td>
      <td>5</td> <td>2</td> <td>1..6</td> <td>1..9</td> <td>-</td> <td>-</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>O</td> <td>O</td>
      </tr>
  
  <tr><td>Reverse biorthogonal spline </td>
      <td>6</td> <td>2</td> <td>1..9</td> <td>1..6</td> <td>-</td> <td>-</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>O</td> <td>O</td> 
      </tr>
  
  <tr><td>Discrete Meyer </td>
      <td>7</td> <td>3</td> <td>-</td><td>-</td> <td>8</td> <td>>=0</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Meyer </td>
      <td>8</td> <td>3</td> <td>-</td> <td>-</td> <td>8</td> <td>>=0</td>  <td>-</td> <td>-</td> <td>O</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Gaussian </td>
      <td>9</td> <td>4</td> <td>1..8</td> <td>-</td><td>5</td> <td>>=0</td>  <td>-</td> <td>-</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Mexican hat </td>
      <td>10</td> <td>4</td> <td>-</td> <td>-</td><td>5</td> <td>>=0</td>  <td>-</td> <td>-</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Morlet </td>
      <td>11</td> <td>4</td> <td>-</td> <td>-</td><td>4</td> <td>>=0</td>  <td>-</td> <td>-</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Complex Gaussian </td>
      <td>12</td> <td>5</td> <td>1..8</td><td>-</td><td>5</td> <td>>=0</td>  <td>-</td> <td>-</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Complex Morlet </td>
      <td>13</td> <td>5</td> <td>-</td> <td>-</td><td>4</td> <td>>=0</td>  <td>>0</td> <td>>0</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Complex Shannon </td>
      <td>14</td> <td>5</td> <td>-</td> <td>-</td><td>20</td> <td>>=0</td>  <td>>0</td> <td>>0</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
  <tr><td>Complex frequency B-spline </td>
      <td>15</td> <td>5</td> <td>1,2,...</td><td>-</td> <td>20</td> <td>>=0</td>  <td>>0</td> <td>>0</td> <td>-</td> <td>O</td> <td>-</td> <td>-</td>
      </tr>
  
</table>

<h4><font color=\"#008000\">Available orders for biorthogonal and reverse biorthogonal wavelets</font></h4>
<p>
For biorthogonal and reverse biorthogonal wavelets, only some combinations of the orders are available. For details, refer to the 
functions, <a href=\"modelica://Wavelet.Families.wavBiorSpline\">wavBiorSpline</a> and <a href=\"modelica://Wavelet.Families.wavRevBiorSpline\">wavRevBiorSpline</a>.
</p>

<h4><font color=\"#008000\">Wavelet types</font></h4>
<p>
All wavelets are categorized into five types. They are defined as follows:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Type</b></td>
      <td><b>Description</b></td>
      </tr>
  
  <tr><td>1</td> <td>Orthogonal wavelets</td> </tr>
  <tr><td>2</td> <td>Biorthogonal wavelets</td> </tr>
  <tr><td>3</td> <td>Non-(bi)orthogonal wavelets with scaling function</td> </tr>
  <tr><td>4</td> <td>Non-(bi)orthogonal wavelets without scaling function</td> </tr>
  <tr><td>5</td> <td>Complex wavelets (without scaling function)</td> </tr>
</table>


</html>"));

  end Families;

  package General "General functions"

    function _fft "An external function to carry out FFT"
        input Integer dir = 1
        "Calculation direction: 1 - forward; -1 - backward";
        input Integer m = 3 "Data points = 2^m";
        input Real[:] x = {i for i in 1:8} "Real part of the data";
        input Real[:] y = {0,0,0,0,0,0,0,0} "Imaginary part of the data";
        output Real[:] x_out = x "Real part of the result";
        output Real[:] y_out = y "Imaginary part of the result";
        external "C" fft_c(dir, m, x_out, y_out);

        annotation (Include="#include <fft_c.c>");
    end _fft;

    function cumSum "Cumulative sum of a vector. Data type is Real."
        input Real i_v[:]={1,2,3,4} "Input vector";
        output Real o_v[size(i_v,1)] "Output vector";
    algorithm

        o_v[1] := i_v[1];
        // cumulative sum, maybe we should check the type of the wavelets here
        for i in 2:size(i_v,1) loop
            o_v[i] :=o_v[i-1] + i_v[i];
        end for;

    end cumSum;

  function cumSumInt "Cumulative sum of a vector. Data type is Integer"
      input Integer i_v[:]={1,2,3,4} "Input vector";
      output Integer o_v[size(i_v,1)] "Output vector";
  algorithm

      o_v[1] := i_v[1];
      // cumulative sum, maybe we should check the type of the wavelets here
      for i in 2:size(i_v,1) loop
          o_v[i] :=o_v[i-1] + i_v[i];
      end for;

  end cumSumInt;

    function diff "Difference between every two adjacent elements of a vector"
      input Real[:] in_v= {1,2,3,4} "Input vector";
      output Real[size(in_v, 1) - 1] out_v "Output vector";

    algorithm
      for i in 1:(size(in_v, 1) - 1) loop
        out_v[i] := in_v[i+1] - in_v[i];
      end for;

    end diff;

    function fft "Fast Fourier transform"
      import Modelica.ComplexMath.*;

        input Complex[:] x = {Complex(re=1, im=0), Complex(re=2, im=0), Complex(re=3, im=0), Complex(re=4, im=0), Complex(re=5, im=0), Complex(re=6, im=0), Complex(re=7, im=0), Complex(re=8, im=0)}
        "The vector to be transformed";
        input Integer n = 8 "Number of data points to be calculated";
        output Complex[size(x,1)] y "The result";
    protected
        Integer q = integer(log10(n)/log10(2)) "number of bits";
        Real[size(x,1)] x_re;
        Real[size(x,1)] x_im;
        Real[size(x,1)] y_re;
        Real[size(x,1)] y_im;
    algorithm
        for i in 1:n loop
            x_re[i] := real(x[i]);
            x_im[i] := imag(x[i]);
        end for;

        (y_re, y_im) := _fft(1, q, x_re, x_im);

        for i in 1:n loop
            y[i].re := y_re[i];
            y[i].im := y_im[i];
        end for;

    end fft;

    function fftShift "Shift zero-frequency component to center of spectrum"
      import Modelica.ComplexMath.*;

        input Complex[:] x "The vector to be shifted";
        output Complex[size(x,1)] y "The shifted vector";
    protected
        Integer m = size(x,1);
        Integer p;
    algorithm
        p := integer(ceil(m/2));
        y := cat(1, x[{i for i in p+1:m}], x[{i for i in 1:p}]);
    end fftShift;

  function findIndex "Find the location of a value in a monotone vector"
    input Real u[:]={1,2,3,3,3,4,5,6,7,8}
        "A monotone vector, in which u[k+1] >= u[k] for all k";
    input Real x=3 "The value to be located in u";
    output Integer id1 "The largest index, id1, so that u[id1] <= x";
    output Integer id2 "The smallest index, id2, so that u[id2] >= x";
    protected
    Integer k;
  algorithm
    id1:=1;
    id2:=size(u,1);

    if x < u[1] then // special cases
      id1 := 0;
      id2 := 0;

    elseif x > u[id2] then // special cases
      id1 := id2+1;
      id2 := id1;

    else  // search
      while id2-id1>1 loop
        k:=integer((id2-id1)/2)+id1;
        if x<= u[k] then
          id2 := k;
        else
          id1 := k;
        end if;
      end while;

    // check if x is exactely equal to one element
    if x == u[id1] then
      id2 := id1;
    elseif x == u[id2] then
      id1 := id2;
    end if;

    end if;

  annotation (Inline=true, Documentation(info="<html>
<p>
This function searches the index of the input vector, u, according to x, so that:
</p><p>
If x is between u[k] and u[k+1], id1 = k and id2 = k+1.
</p><p>
If x is equal to u[k], id1 = id2 = k.
</p><p>
If x is smaller than all elements of u, id1 = id2 = 0.
</p><p>
If x is larger than all elements of u, id1 = id2 = size(u,1)+1.
</p><p>
The searching method is binary search, which has high efficiency for large searching spaces.
</p>
</html>"));
  end findIndex;

       function filterBank
      "Get the four wavelet filters based on a given scaling filter"
             input Real F[:]={0.5, 0.5} "The scaling filter";
             output Real lod[:] "High pass filter for decomposition";
             output Real hid[:] "Low pass filter for decomposition";
             output Real lor[:] "High pass filter for reconstruction";
             output Real hir[:] "Low pass filter for reconstruction";

    protected
             Real F1[:];

       algorithm
             // normalized filter
             F1 := F/sum(F);

             // generate the filters
             lor :=sqrt(2).*F1;
             hir :=quadReverse(lor);
             lod :=Modelica.Math.Vectors.reverse(lor);
             hid :=Modelica.Math.Vectors.reverse(hir);
       end filterBank;

    function ifft "Inverse fast Fourier transform"
      import Modelica.ComplexMath.*;

        input Complex[:] x = {Complex(re=36, im=0), Complex(re=-4, im=9.65685424949238), Complex(re=-4, im=4), Complex(re=-4, im=1.65685424949238), Complex(re=-4, im=0), Complex(re=-4, im=-1.65685424949238), Complex(re=-4, im=-4), Complex(re=-4, im=-9.65685424949238)}
        "The vector to be inverse transformed";
        input Integer n = 8 "Number of data points for calculation";
        output Complex[size(x,1)] y "The result";
    protected
        Integer q = integer(log10(n) / log10(2)) "number of bits";
        Real[size(x,1)] x_re;
        Real[size(x,1)] x_im;
        Real[size(x,1)] y_re;
        Real[size(x,1)] y_im;
    algorithm

        for i in 1:n loop
            x_re[i] := real(x[i]);
            x_im[i] := imag(x[i]);
        end for;

        (y_re, y_im) := _fft(-1, q, x_re, x_im);

        for i in 1:n loop
            y[i].re := y_re[i];
            y[i].im := y_im[i];
        end for;

    end ifft;

       function innerProduct "Inner product of two same length vectors"
             input Real x1[:]={1,2,3,4} "Vector 1";
             input Real x2[size(x1,1)]={1,2,3,4} "Vector 2";
             output Real y=0. "Result";

    protected
             Integer n;  // vector length
             Integer k;  // loop counter

       algorithm
             n :=size(x1,1);
             for k in 1:n loop
               y :=y + x1[k]*x2[k];
             end for;

       end innerProduct;

       function interpL "One-dimensional linear interpolation"
             input Real u[:]={1,2,3,4} "The variable of the original data";
             input Real y[size(u,1)]={0,1,-1,2}
        "The function of the original data";
             input Real u1[:]={0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5}
        "The variable of the interpolated function";
             output Real y1[size(u1,1)]
        "The interpolated result at the variable values, u1[:]";

    protected
             Integer nu "length of u";
             Integer nu1 "length of u1";
             Integer ku "loop counter for nu";
             Integer ku1 "loop counter for nu1";

       algorithm
             // get the vector length
             nu := size(u,1);
             nu1 := size(u1,1);

             // interpolation
             ku := 1;
             for ku1 in 1:nu1 loop
               // u1-elements left to u-range
               if u1[ku1] < u[1] then
                 y1[ku1] := y[1];

               // u1-elements in the u-range
               elseif u1[ku1] < u[nu] then
                 while ku < nu loop
                   if u1[ku1] >= u[ku] and u1[ku1] < u[ku+1] then
                     y1[ku1] := y[ku] + (y[ku+1] - y[ku]) * (u1[ku1] - u[ku]) / (u[ku+1] - u[ku]);
                     break;
                   else
                     ku := ku + 1;
                   end if;
                 end while;

               // u1-elements right to u-range
               else
                 y1[ku1] :=y[nu];
               end if;
             end for;

             annotation (Inline=true, Documentation(info="<html>
<p>
This function finds the values at u1 of the linearly interpolated function of the original data pairs, u and y. 
Both u and u1 must be monotone, meaning, e.g. u[k]>=u[k-1]. For the elements in u1 whose values are less than u[1], 
the corresponding elements in y1 are set as y[1]. For the elements in u1 whose values are greater than the last 
value of u, the corresponding elements in y1 are set as y[end]. 
</p>
</html>"));
       end interpL;

    function midVector "Extract the middle part of a vector"
      input Real x[:]={1,2,3,4,5,6} "Original vector";
      input Integer n=3 "Number of the vector elements to be extracted";
      output Real y[:] "The result vector";

    protected
      Integer nx=size(x, 1) "length of vector x";
      Integer id1 "first element number";
      Integer id2 "last  element number";

    algorithm
      if n >= nx then   // check the input
          y := x;

      elseif mod(nx,2) == mod(n, 2) then
          id1 := integer((nx - n)/2 +1);
          id2 := nx - id1 +1;
          y := x[id1:id2];

      else
          id1 := integer(floor((nx - n)/2) +1);
          id2 := nx - id1;
          y := x[id1:id2];

      end if;

    end midVector;

    function nStdIfft "Inverse non-standard 1-D fast Fourier transform"
      import Modelica.ComplexMath.*;

        //input Real[:] xhat = {0, 0, 0, 0, 0, 0, 0.00979994304697744, 0.222461474833219, 0.707106781186548, 0.974941481431080, 0.999951979405149, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999951979405149, 0.974941481431080, 0.707106781186548, 0.222461474833219, 0.00979994304697744, 0, 0, 0, 0, 0}
        input Complex[:] xhat "The vector to be inverse transformed";
        input Real low = -8 "Lower boundary of the function";
        input Real high = 8 "Higher boundary of the function";
        output Real[size(xhat,1)] x "The recovered signal";
        output Real[size(xhat,1)] t "Time grid";
    protected
        Real pi = Modelica.Constants.pi;
        Integer n = size(xhat,1);
        Real delta;
        Real[size(xhat,1)] omega;
        Complex[size(xhat,1)] tmp1;
        Complex[size(xhat,1)] tmp2;
        Complex[size(xhat,1)] x_c;
    algorithm
        // time grid resolution
        delta := (high - low)/n;
        // frequency grid
        omega := {i for i in (-n):2:(n-2)} / (2*n*delta);
        tmp1 := exp(2*pi*j*low*omega);
        for i in 1:n loop
            tmp2[i] := xhat[i] * tmp1[i] / delta;
        end for;

        tmp2 := fftShift(tmp2);

        x_c := ifft(tmp2, n);

        for i in 1:n loop
            // maybe we have to check the imaginary part here?
            x[i] := x_c[i].re;
        end for;

        t := low .+ {i for i in 0:n-1}.*delta;

    end nStdIfft;

       function quadReverse "Quadrature mirror of a given vector"

             input Real x[:]={1,2,3,4,5,6,7,8} "Vector";
             input Boolean p = true
        "If p is true, the signs of the even index entries of the reversed input vector, x, is reversed. Otherwise, odd index entries are reversed.";
             output Real y[size(x, 1)]
        "Quadrature mirrored elements of the input vector";

    protected
             Integer id1 = 0;
             Integer id2 = 0;

       algorithm
             y := Modelica.Math.Vectors.reverse(x);

             // start index number to inverse the sign
             if p then
               id1:=2;
             else
               id1:=1;
             end if;

             // half of the elements sign-inversed
             for id2 in id1:2:size(y,1) loop
               y[id2] :=-1*y[id2];
             end for;

       annotation (Inline=true, Documentation(info="<html>
<p>
Changes the signs of the even index entries of the reversed input vector, x, if p is true. 
If p is false the same holds for odd index entries.
</p>
</html>"));
       end quadReverse;

    function sinc "Sinc function"
      import Modelica.Constants.*;
      input Real x=pi/2;
      output Real y;

    algorithm
        if x==0 then
          y:=1;
        else
          y :=sin(pi*x)/(pi*x);
        end if;

    end sinc;

       function upsample
      "Up-sampling of a vector (insert a zero after every element except the last one)"
             input Real u[:]={1,2,3,4} "Input vector";
             output Real y[size(u,1)*2-1]
        "Output vector with the (double length of u) -1";

    protected
             Integer nu "input vector length";
             Integer k "counter";

       algorithm
             nu :=size(u, 1);
             y :=fill(0., 2*nu-1);
             for k in 1:nu loop
               y[k*2-1] :=u[k];
             end for;

       end upsample;

       function wavConv
      "Fully convolving of a data vector and a filter vector for wavelet transform"
             input Real u[:]={1,2,3,4} "The data vector";
             input Real f[:]={1,2} "The filter vector";
             output Real y[size(u,1)+size(f,1)-1] "The filtered data vector";

    protected
             Real u1[size(u,1)+2*(size(f,1)-1)] "the extended vector";
             Real f1[size(f,1)] "the reversed filter";
             Integer nu1 "length of the extended vector x1";
             Integer nf "filter length";
             Integer k "loop counter";

       algorithm
             // vector lengths
             nu1 :=size(u,1)+2*(size(f,1)-1);
             nf :=size(f, 1);

             // extend the input data vector using zero padding
             u1 :=cat(
             1,
             fill(0, nf - 1),
             u,
             fill(0, nf - 1));

             // the mirrored filter
             f1 :=Modelica.Math.Vectors.reverse(f);

             // convolution
             for k in 1:nu1-nf+1 loop
               y[k] :=innerProduct(u1[k:k + nf - 1], f1);
             end for;

       end wavConv;

  annotation (Documentation(info="<html>
<p>
This section defines some general functions, which are common utilities for wavelet transformation and other functions.
</p>
</html>"));
  end General;

  package Records "Records, mainly used for graphic user interface"

    record denoisParameters "Parameters for 1D wavelet denoising"
      extends Modelica.Icons.Record;

        Integer decLevel=2 "Decomposition levels, <=12" annotation(Dialog(group="Parameters"));
        Types.threshMethod thMethod=1
        "Select what method is to be used for calculating the threshold" annotation(Dialog(group="Parameters"));
        Boolean sorh=true "Soft (false) or hard (true) thresholding method"
                                                                           annotation(Dialog(group="Parameters"));

        Real thD12=-1 "Detail level 12" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD11=-1 "Detail level 11" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD10=-1 "Detail level 10" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD9=-1 "Detail level 9" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD8=-1 "Detail level 8" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD7=-1 "Detail level 7" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD6=-1 "Detail level 6" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD5=-1 "Detail level 5" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD4=-1 "Detail level 4" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD3=-1 "Detail level 3" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD2=-1 "Detail level 2" annotation(Dialog(group="Thresholds (>=0) for all levels"));
        Real thD1=-1 "Detail level 1" annotation(Dialog(group="Thresholds (>=0) for all levels"));

    annotation (Documentation(info="<html>
<p>
Parameters for carrying out wavelet denoising. 
</p>

</html>"));

    end denoisParameters;

    record mraParameters "Parameters for multi-resolution analysis (MRA)"
    extends Modelica.Icons.Record;

        Integer decLevel=2 "Decomposition levels, <=12" annotation(Dialog(group="Wavelet"));
        Types.mraDisplay mraData=2
        "Select what data type is to be displayed in the decomposed levels"                             annotation(Dialog(group="Wavelet"));

        Real rA=1 "Approximating coefficients" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD12=1 "Detail coefficients, level 12" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD11=1 "Detail coefficients, level 11" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD10=1 "Detail coefficients, level 10" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD9=1 "Detail coefficients, level 9" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD8=1 "Detail coefficients, level 8" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD7=1 "Detail coefficients, level 7" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD6=1 "Detail coefficients, level 6" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD5=1 "Detail coefficients, level 5" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD4=1 "Detail coefficients, level 4" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD3=1 "Detail coefficients, level 3" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD2=1 "Detail coefficients, level 2" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));
        Real rD1=1 "Detail coefficients, level 1" annotation(Dialog(group="Tuning factor of the coefficients before reconstruction. Levels higher than decLevel are omitted."));

    annotation (Documentation(info="<html>
<p>
Parameters for carrying out wavelet multi-resolution analysis (MRA). 
</p>

</html>"));

    end mraParameters;

    record wavFuncOut "Output data of function, wavFunc()"
        String wavName "Wavelet name" annotation(Dialog(group="Output"));
        Integer wavType "Wavelet type: 1. orthogonal wavelets; 2. biorthogonal wavelets; 
    3. wavelet with scaling function; 4. wavelet without scaling function; 5. complex wavelet without scaling function"
                                                                                                            annotation(Dialog(group="Output"));
        Real F1[:]
        "Orthogonal wavelets: scaling filter; biorthogonal wavelets: scaling filter for decomposition"
                                                                                                      annotation(Dialog(group="Output"));
        Real F2[:]
        "Only for biorthogonal wavelets: scaling filter for reconstruction"
                                                                           annotation(Dialog(group="Output"));
        Real lod[:] "Low pass filter for decomposition"
                                                       annotation(Dialog(group="Output"));
        Real hid[:] "High pass filter for decomposition"
                                                        annotation(Dialog(group="Output"));
        Real lor[:] "Low pass filter for reconstruction"
                                                        annotation(Dialog(group="Output"));
        Real hir[:] "High pass filter for reconstruction"
                                                         annotation(Dialog(group="Output"));
        Real x[:] "Equidistance grid"
                                     annotation(Dialog(group="Output"));
        Real phi1[:] "Scaling functions"
                                        annotation(Dialog(group="Output"));
        Real psi1[:] "Wavelet functions"
                                        annotation(Dialog(group="Output"));
        Real phi2[:]
        "Scaling functions, only applicable for biorthogonal wavelets"             annotation(Dialog(group="Output"));
        Real psi2[:]
        "Wavelet functions, only applicable for biorthogonal wavelets"             annotation(Dialog(group="Output"));
    annotation (Documentation(info="<html>
<p>
Output parameters of function <a href=\"modelica://Wavelet.Families.wavFunc\">wavFunc()</a>.
</p>

</html>"));
    end wavFuncOut;

    record wavletDefinition "Wavelet definition"
    extends Modelica.Icons.Record;
        Types.wavletID wavID=2
        "ID number of the wavelet family to be displayed, 1..15."                        annotation(Dialog(group="Input"));
        Integer Nd=3 "Wavelet order" annotation(Dialog(group="Input"));
        Integer Nr=3
        "Wavelet order for reconstruction, only used for biorthogonal wavelets"
                                                                                               annotation(Dialog(group="Input"));
        Real fb=1.0 "Bandwidth frequency, only for complex wavelets" annotation(Dialog(group="Input"));
        Real fc=1.0 "Center frequency, only for complex wavelets" annotation(Dialog(group="Input"));
        Real range=5.0
        "Boundary of the wavelet and scaling functions, [-range, range]. Valid for wavelets with theoretically unlimited support. Effective support is used as default value if range==0"
                                                                                                            annotation(Dialog(group="Input"));
        Integer refinement=6 "Refinement of the estimated wavelet and scaling functions. 
    For orthogonal and biorthogonal wavelets, it stands for the iteration 
    times to generating functions. The returned function length is (filter_length-1)*2^refinement+1. 
    For other wavelets, the length of the returned functions is 2^refinement."     annotation(Dialog(group="Input"));

    annotation (Documentation(info="<html>
<p>
Parameters for wavelet definition. 
</p>
<p>
Refer to the section of wavelet <a href=\"modelica://Wavelet.Families\">Families</a> 
for detailed information about the available wavelets.

</p>

</html>"));
    end wavletDefinition;

  end Records;

  package Types "Definitions of enumeration values"

    type mraDisplay = enumeration(
        coefficients "1, Wavelet coefficients",
        signal "2, Reconstructed signal")
      "Selection of display data for MRA: 1 - wavelet coefficients; 2 - reconstructed signal";

  type threshMethod = enumeration(
        fixedForm "1_fixedForm",
        SURE "2_SURE",
        heurSURE "3_heurSURE",
        miniMax "4_MiniMax") "Methods for calculating threshold for denoising: 1 - fixed form; 2 - SURE method; 
  3 - heuristic SURE; 4 - minimal maximum method";
    type wavletID = enumeration(
        Haar "1_Haar",
        Daubechies "2_Daubechies",
        Symlets "3_Symlets",
        Coiflets "4_Coiflets",
        BiorSpline "5_BiorSpline",
        RBiorSpline "6_RBiorSpline",
        DiscreteMeyer "7_DiscreteMeyer",
        Meyer "8_Meyer",
        Gaussian "9_Gaussian",
        MexicanHat "10_MexicanHat",
        Morlet "11_Morlet",
        XGaussian "12_XGaussian",
        XMorlet "13_XMorlet",
        XShannon "14_XShannon",
        XFBSpline "15_XFBSpline") "Definition of wavelet identifiers: 1 - Haar; 2 - Daubechies; 3 - Symlets; 4 - Coiflets; 
  5 - Biorthogonal spline; 6 - Reverse biorthogonal spline; 7 - Meyer; 8 - Discrete Meyer; 9 - Gaussian; 10 - Mexican hat; 
  11 - Morlet; 12 - Complex Gaussian; 13 - Complex Morlet; 14 - Complex Shannon; 15 - Complex frequency B-Spline";

  annotation (Documentation(info="<html>
<p>Enumeration variables are defined here. 

<h4><font color=\"#008000\">waveletID: ID numbers of wavelets</font></h4>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>waveletID</b></td>
      <td><b>Wavelet name</b></td>
      </tr>
  
  <tr><td>1</td> <td>Haar</td> </tr>
  <tr><td>2</td> <td>Daubechies</td> </tr>
  <tr><td>3</td> <td>Symlets</td> </tr>
  <tr><td>4</td> <td>Coiflets</td> </tr>
  <tr><td>5</td> <td>Biorthogonal spline</td> </tr>
  <tr><td>6</td> <td>Reverse biorthogonal spline</td> </tr>
  <tr><td>7</td> <td>Meyer</td> </tr>
  <tr><td>8</td> <td>Discrete Meyer</td> </tr>
  <tr><td>9</td> <td>Gaussian</td> </tr>
  <tr><td>10</td> <td>Mexican hat</td> </tr>
  <tr><td>11</td> <td>Morlet</td> </tr>
  <tr><td>12</td> <td>Complex Gaussian</td> </tr>
  <tr><td>13</td> <td>Complex Morlet</td> </tr>
  <tr><td>14</td> <td>Complex Shannon</td> </tr>
  <tr><td>15</td> <td>Complex frequency B-Spline</td> </tr>
  
</table>

<h4><font color=\"#008000\">mraDisplay: Data type to be displayed in MRA</font></h4>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Value</b></td>
      <td><b>Name</b></td>
      <td><b>Description</b></td>
      </tr>
  
  <tr><td>1</td> <td>coefficients</td> <td>Wavelet coefficients</td> </tr>
  <tr><td>2</td> <td>signal</td> <td>Reconstructed signal</td> </tr>
  
</table>

<h4><font color=\"#008000\">threshID: ID numbers of different methods for threshold calculation for denoising</font></h4>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Value</b></td>
      <td><b>Name</b></td>
      <td><b>Description</b></td>
      </tr>
  
  <tr><td>1</td> <td>fixedForm</td> <td>Fixed form threshold. Threshold value = sqrt(2*ln(N)), with N being the length of the data to be denoised.</td> </tr>
  <tr><td>2</td> <td>SURE</td> <td>The threshold is calculated using the principle of Stein's Unbiased Risk Estimate(SURE). </td> </tr>
  <tr><td>3</td> <td>heurSure</td> <td>Heuristic SURE method. It is a combination of fixedForm and SURE methods. The noise level is firstly tested. 
  For a high signal-to-noise ratio, SURE method is used, otherwise fixedForm method is used. </td> </tr>
  <tr><td>4</td> <td>miniMax</td> <td>The threshold is calculated according to the miniMax principle, which realizes the minimum of the maximum mean square error.</td> </tr>
    
</table>
    
</p>
</html>"));
  end Types;

annotation (Documentation(info="<html>
<p>
The package <b>Wavelet</b> is a function library developed for Modelica to carry out wavelet transform and related calculations.
It is used for the post processing of the simulation results or other data. Online operation, meaning execution of the functions 
during the simulation process, is not possible.
</p><p>
The numerical calculation of wavelet transform requires that the data be sampled in equidistant time grids. 
If non-equidistant time grids are used, the data should be firstly converted to equidistant using function 
<a href=\"modelica://Wavelet.General.interpL\">interpL</a> in this library. If the data are to be generated by Dymola 
with model simulation, it is suggested that <b>Equidistant time grid</b> should be checked in the <b>Simulation</b> window through the menu
 <b>Simulation/Setup.../Output/Output selection</b>.
</p><p>
Most functions in this library do not check the input parameters. It is the task of the user to provide valid parameters. 
Upon wrong input parameters, unexpected errors might occur.
</p><p>
This library is developed and tested with the demo version of Dymola 2013.
</p><p>
<h4><font color=\"#008000\">List of the packages:</font></h4>
<p>Detailed information is to be found in each package.</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Package</b></td>
      <td><b>Description</b></td>
      </tr>
  
  <tr><td><a href=\"modelica://Wavelet.Examples\">Examples</a></td> <td>Some examples showing some functionalities of this library</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.MRA\">MRA</a></td> <td>Wavelet application for multi-resolution analysis (MRA)</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.Denoising\">Denoising</a></td> <td>Wavelet application for denoising</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.Transform\">Transform</a></td> <td>Wavelet transform and tightly related functions</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.Families\">Families</a></td> <td>Definition of wavelet families</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.General\">General</a></td> <td>General purpose functions used in this library</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.Records\">Records</a></td> <td>Definitions for supporting graphic user interface</td> </tr>
  <tr><td><a href=\"modelica://Wavelet.Types\">Types</a></td><td>Definitions of enumeration values </td> </tr>
</table>
</p>

<br>
<h4><font color=\"#008000\">Required libraries:</font></h4>
Besides the Modelica standard conform library, the following libraries are required for using this wavelet library:
</p>
<p style=\"margin-left: 40px;\">
<b>* Modelica_LinearSystems2:</b> This library is used for displaying the curves.<br>
<b>* Plot3D:</b> This is used for showing the images generated by wavelet continuous transform (CWT). 
The library Plot3D is delivered with Dymola. In the Dymola demo version, the functionality of Plot3D 
seems to be limited, such that the 3D surfaces can only be displayed in a 2D manner.
<br><p>

<h4><font color=\"#008000\">License:</font></h4>
<p>All files in this directory (Wavelet) and in all subdirectories, especially all files that build package \"Wavelet\" 
are released under the <b>Modelica License 2</b>.</p>
</p>
<br>
<h4><font color=\"#008000\">Developers:</font></h4>
<p style=\"margin-left: 40px;\">
Dr. Michael Gao, michael.gao@tum.de<br>
Mr. Qipeng Hu, qipeng.hu@gmail.com<br>
Mr. Weihua Wang<br>
Mr. Hanko Ipach</p>
<h4><font color=\"#008000\">Acknowlegement:</font></h4>
The development of this library is supported by EU within the project, Clean Sky, sub-project, MoMoLib, No. 296369.
The first delivery of this Wavelet library is made in October 2013.


</p>
<br>
</html>"),    uses(Modelica(version="3.2"), Complex(version="1.0")));
end Wavelet;
