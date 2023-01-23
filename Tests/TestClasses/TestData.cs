using System.Runtime.CompilerServices;
using Chemistry;
using MathNet.Numerics.Distributions;
using Microsoft.VisualStudio.TestPlatform.CommunicationUtilities;
using Plotly.NET;
using ThermoFisher.CommonCore.Data.Business;
using static Microsoft.FSharp.Core.ByRefKinds;

namespace Tests;
public abstract class TestData
{
    public double[] Xarray { get; set; }
    public double[] Yarray { get; set; }
    public int Length { get; set; }
    private double Spacing { get; set; }
    protected TestData(int length, double startValue, double spacing)
    {
        Length = length;
        Spacing = spacing;

        Xarray = Enumerable.Range(0, length)
            .Select(i => i * spacing + startValue)
            .ToArray();
        Yarray = Enumerable.Range(0, length)
            .Select(i => i * spacing + startValue)
            .ToArray();
    }

    protected TestData()
    {

    }

    /// <summary>
    /// Applies a function to each element in an array. 
    /// </summary>
    /// <param name="function"></param>
    public void ApplyElementwise(Func<double, double> function)
    {
        for (int i = 0; i < Length; i++)
        {
            Yarray[i] = function(Yarray[i]);
        }
    }

    public GenericChart.GenericChart Plot()
    {
        return Chart2D.Chart.Scatter<double, double, int>(Xarray, Yarray, StyleParam.Mode.Lines);
    }

    public static TestData operator +(TestData a, TestData b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
        if (Math.Abs(a.Spacing - b.Spacing) > 0.001) throw new ArgumentException("Incompatible spacing in data");

        for (int i = 0; i < a.Length; i++)
        {
            a.Yarray[i] += b.Yarray[i];
        }

        return a;
    }

    public void ElementwiseArrayAddition(TestData peak)
    {
        for (int i = 0; i < peak.Length; i++)
        {
            Yarray[i] += peak.Yarray[i]; 
        }
    }
    public static TestData operator *(TestData a, TestData b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
        if (Math.Abs(a.Spacing - b.Spacing) > 0.001) throw new ArgumentException("Incompatible spacing in data");

        for (int i = 0; i < a.Length; i++)
        {
            a.Yarray[i] *= b.Yarray[i];
        }

        return a;
    }
    public static TestData operator -(TestData a, TestData b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
        if (Math.Abs(a.Spacing - b.Spacing) > 0.001) throw new ArgumentException("Incompatible spacing in data");

        for (int i = 0; i < a.Length; i++)
        {
            a.Yarray[i] -= b.Yarray[i];
        }

        return a;
    }
    public static TestData operator *(TestData a, double scalar)
    {
        for (int i = 0; i < a.Length; i++)
        {
            a.Yarray[i] *= scalar;
        }

        return a;
    }
}

public class TwoPeakTestData
{
    public ChromatographicPeak Chromatogram { get; set; }
    public List<GaussianPeak> Data { get; set; }
    private Normal ShotNoiseDistribution { get; set; }
    public int NumberSpectra { get; }
    public (double, double) Peak1 { get; }
    public (double, double) Peak2 { get; }
    public List<(double Peak1, double Peak2, double ChromRatio)> PeakRatios { get; set; }
    public int Length { get; }
    public TwoPeakTestData(int numberSpectra, int length, (double, double) peak1,
        (double, double) peak2, ChromatographicPeak chromatogram,
        Normal shotNoiseDistribution)
    {
        Data = new();
        Chromatogram = chromatogram;
        ShotNoiseDistribution = shotNoiseDistribution;
        NumberSpectra = numberSpectra;
        Peak1 = peak1;
        Peak2 = peak2;
        Length = length;
        PeakRatios = new();
    }
    public void CreateDataSet(Normal noiseDistribution, double peak1BaseIntensity,
        double peak2BaseIntensity)
    {
        for (int i = 0; i < NumberSpectra; i++)
        {
            // generate peak1 and peak2
            // multiply peak1 and peak2 by their chromatogram multiplier
            // multiply by the shot noise distribution. 
            // add the spectra together
            // add to list 

            double peak1ShotMultiple = ShotNoiseDistribution.Sample();
            double peak2ShotMultiple = ShotNoiseDistribution.Sample();

            double chromRatio = Chromatogram.Yarray[i] / Chromatogram.Max;

            double intensityPeak1 = chromRatio * (peak1BaseIntensity + peak1ShotMultiple * peak1BaseIntensity);
            double intensityPeak2 = chromRatio * (peak2BaseIntensity + peak2ShotMultiple * peak2BaseIntensity);


            PeakRatios.Add((intensityPeak1, intensityPeak2, chromRatio));

            GaussianPeak peak1 = new(Peak1.Item1, Peak1.Item2,
                Length, 0, 1.0, noiseDistribution, intensityPeak1);
            GaussianPeak peak2 = new(Peak2.Item1, Peak2.Item2,
                Length, 0, 1.0, null, intensityPeak2);

            Data.Add(peak2 + peak1 as GaussianPeak);
        }
    }

    public (double[][] xArrays, double[][] yArrays) GetJaggedArrays()
    {
        double[][] xArrays = new double[Data.Count][];
        double[][] yArrays = new double[Data.Count][];

        for (int i = 0; i < Data.Count; i++)
        {
            xArrays[i] = Data[i].Xarray;
            yArrays[i] = Data[i].Yarray;
        }
        return (xArrays, yArrays);
    }

    public double[] CalculateTics()
    {
        return Data.Select(i => i.Yarray.Sum()).ToArray();
    }
}

public class SimulatedProtein : GaussianPeak
{
    IsotopicDistribution _distribution;


    /// <summary>
    /// Record type used to facilitate bin, mz, and intensity matching. 
    /// </summary>
    /// <param name="Bin">Integer bin number.</param>
    /// <param name="Mz">Mz value.</param>
    /// <param name="Intensity">Intensity value.</param>
    internal readonly record struct BinValue(int Bin, double Mz, double Intensity);

    /// <summary>
    /// Custom comparer to use in override for List.BinarySearch() that accepts a custom comparer as an argument. 
    /// </summary>
    internal class BinValueComparer : IComparer<BinValue>
    {
        public int Compare(BinValue x, BinValue y)
        {
            return x.Bin.CompareTo(y.Bin);
        }
    }
    private static List<BinValue> CreateBinValueList(double[] xArray, double[] yArray,
        double min, double binSize)
    {
        var binIndices = xArray
            .Select((w, i) =>
                new { Index = i, Bin = (int)Math.Round((w - min) / binSize, MidpointRounding.AwayFromZero) });
        List<BinValue> binValues = new List<BinValue>();
        foreach (var bin in binIndices)
        {
            binValues.Add(new BinValue(bin.Bin, xArray[bin.Index], yArray[bin.Index]));
        }

        return binValues;
    }

    private (List<double> xVals, List<double> yVals) ConsumeDistribution(IsotopicDistribution distribution, double binSize)
    {
        double min = distribution.Masses.Min() - 1; 

        int numberOfBins = Length;
        List<BinValue> listBinValues = CreateBinValueList(distribution.Masses.ToArray(),
            distribution.Intensities.ToArray(), min, binSize);
        listBinValues.Sort((m, n) => m.Bin.CompareTo(n.Bin));
        
        List<double> xVals = new();
        List<double> yVals = new();
        
        for (int i = 0; i < numberOfBins; i++)
        {
            List<BinValue> binValRecord = new();
            int index = listBinValues.BinarySearch(new BinValue() { Bin = i }, new BinValueComparer());
            if (index < 0)
            {
                index = ~index;
            }

            int k = 0;
            while (k + index <= listBinValues.Count - 1)
            {
                // binary search gets just the first index that it finds, so you 
                // need to check on the left and right side for elements that 
                // match the bin value. 

                if (index + k < listBinValues.Count && listBinValues[index + k].Bin == i)
                {
                    binValRecord.Add(listBinValues[index + k]);
                    if (k == 0)
                    {
                        k++;
                        continue;
                    }
                }

                if (index - k > 0 && listBinValues[index - k].Bin == i)
                {
                    binValRecord.Add(listBinValues[index - k]);
                }

                k++;
            }

            if (binValRecord.Count() > 1)
            {
                xVals.Add(binValRecord.Average(m => m.Mz));
                yVals.Add(binValRecord.Average(m => m.Intensity));
                continue;
            }

            if (binValRecord.Count == 0)
            {
                // TODO: x array values aren't handled correctly. Should be the value of the m/z axis at the particular bin, not a zero. 
                xVals.Add(0);
                yVals.Add(0);
                continue;
            }

            xVals.Add(binValRecord.First().Mz);
            yVals.Add(binValRecord.First().Intensity);
        }

        return (xVals, yVals); 
    }
    public SimulatedProtein(IsotopicDistribution distribution, int length, double spacing) : base()
    {
        Length = length; 
        _distribution = distribution;
        this.Xarray = new double[length];
        this.Yarray = new double[length];
        Xarray[0] = distribution.Masses.Min();
        for (int i = 1; i < Xarray.Length; i++)
        {
            this.Xarray[i] = spacing * i + Xarray[0];
        }
        // put the isotopic distribution into the x array and the y array 
        var arrays = ConsumeDistribution(distribution, spacing);

        base.Yarray = arrays.yVals.ToArray();
    }
}


public class GaussianPeak : TestData
{
    public double Mean { get; set; }
    public double Stddev { get; set; }
    private Normal PeakDistribution { get; set; }
    public Normal? NoiseDistribution { get; set; }
    public double IntensityMultiple { get; }
    public GaussianPeak(double mean, double stddev,
        int length, double startValue, double spacing,
        Normal? noiseDistribution, double intensity = 1) : base(length, startValue, spacing)
    {
        Mean = mean;
        Stddev = stddev;
        PeakDistribution = new Normal(mean, stddev);
        IntensityMultiple = intensity;

        ApplyElementwise(GaussianFunc);
        if (noiseDistribution != null)
        {
            NoiseDistribution = noiseDistribution;
            AddNoise();
        }
    }

    protected GaussianPeak()
    {

    }

    public void AddNoise(Normal noiseDistribution)
    {
        NoiseDistribution = noiseDistribution;
        AddNoise();
    }
    public virtual void AddShotPeaks(int lowLimitNumberPeaks, int highLimitNumberPeaks, 
        double peakLocationLow, double peakLocationHigh, double peakWidthLow, double peakWidthHigh, 
        double peakIntensityLow, double peakIntensityHigh, int length = 1000, double startValue = 0.0, 
        double spacing = 1.0)
    {
        Random random = new Random();

        // generate a random number of peaks 
        int numberOfRandomPeaks = random.Next(lowLimitNumberPeaks, highLimitNumberPeaks);
        for (int j = 0; j < numberOfRandomPeaks; j++)
        {
            Random randomInnerLoop = new();

            double peakLocation = NextDouble(randomInnerLoop, peakLocationLow, peakLocationHigh);
            // re-roll the peak location if the rolled location is in the mean plusorminus the stddev 
            while (peakLocation > (Mean - Stddev) && peakLocation < (Mean + Stddev))
            {
                peakLocation = NextDouble(randomInnerLoop, peakLocationLow, peakLocationHigh); 
            }
            
            double peakWidth = NextDouble(randomInnerLoop, peakWidthLow, peakWidthHigh);
            double peakIntensity = NextDouble(randomInnerLoop, peakIntensityLow, peakIntensityHigh);
            GaussianPeak newPeak = new GaussianPeak(peakLocation, peakWidth, 
                length, startValue, spacing, null, peakIntensity); 
            ElementwiseArrayAddition(newPeak);
        }
    }

    private double NextDouble(Random random, int low, int high)
    {
        return random.NextDouble() * ((double)high - (double)low) + (double)low;
    }

    private double NextDouble(Random random, double low, double high)
    {
        return random.NextDouble() * (high - low) + low; 
    }
    private void AddNoise()
    {

        ApplyElementwise(AddNoiseFunc);
    }

    protected double AddNoiseFunc(double d) => d + NoiseDistribution.Sample();
    protected double GaussianFunc(double d) => IntensityMultiple * Normal.PDF(Mean, Stddev, d);
}
/// <summary>
/// Derived from Gaussian peak and forced to having a spacing of 1.0d.
/// </summary>
public class ChromatographicPeak : GaussianPeak
{
    public double Max => Yarray.Max();
    public ChromatographicPeak(double timeOfMax, double peakWidth,
        int length, double startValue, Normal? noiseDistribution = null, double spacing = 1d)
        : base(timeOfMax, peakWidth, length, startValue, spacing, noiseDistribution)
    {

    }

}

// generate a new noise data in each iteration of the loop. 
public class NoiseData : TestData
{
    public NoiseData(int length, double startValue, double spacing, double mean, double stddev) : base(length,
        startValue, spacing)
    {
        Normal distribution = new Normal(mean, stddev);
        for (int i = 0; i < Yarray.Length; i++)
        {
            Yarray[i] = distribution.Sample();
        }
    }
    public NoiseData(int length, double startValue, double spacing, double mean, double stddev, Random randomGenerator) : base(length,
        startValue, spacing)
    {
        Normal distribution = new Normal(mean, stddev, randomGenerator);
        for (int i = 0; i < Yarray.Length; i++)
        {
            Yarray[i] = distribution.Sample();
        }
    }

    public void AddNoiseToXArray(Normal distribution)
    {
        for (int i = 0; i < Xarray.Length; i++)
        {
            Xarray[i] += distribution.Sample();
        }
    }
}

