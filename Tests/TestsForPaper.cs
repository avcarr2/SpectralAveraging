using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.Statistics;
using Microsoft.VisualStudio.TestPlatform.CommunicationUtilities;
using NUnit.Framework.Interfaces;
using Plotly.NET;
using Plotly.NET.LayoutObjects;
using SpectralAveraging;
using SpectralAveraging.NoiseEstimates;

namespace Tests
{
    internal class TestsForPaper
    {
        [Test]
        public void TestGeneratePeak()
        {
            Normal noiseDistribution = new Normal(1, 0.01);
            GaussianPeak gaussian = new GaussianPeak(400, 10, 1000, 
                0, 1.0, noiseDistribution); 
            GaussianPeak gaussian2 = new GaussianPeak(600, 10, 1000, 
                0, 1.0, noiseDistribution);
            GaussianPeak addedGauss = (GaussianPeak)(gaussian + gaussian2);
            var plot = addedGauss.Plot(); 
            plot.Show();
        }

        [Test]
        public void TestGenerateChromatographicPeak()
        {
            ChromatographicPeak chroma = new(10, 5, 20, 0, null);
            chroma.Plot().Show(); 
        }

        [Test, Repeat(5)]
        public void TestTwoPeakTestData()
        {
            ChromatographicPeak peak = new(30, 10, 60,
                0, null);
            Normal shotNoiseDistr = new(1, 0.0001);
            Normal thermal = new(8000, 2000); 
            TwoPeakTestData testData = new(60, 1000,
                (400, 10), (600, 10), peak, shotNoiseDistr);
            testData.CreateDataSet(thermal, 100000, 200000);

            var data = testData.GetJaggedArrays();
            var tics = testData.CalculateTics(); 
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.PerformNormalization = true;
            options.BinSize = 1.0;
            options.RejectionType = RejectionType.SigmaClipping;
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.MinSigmaValue = 3.5;
            options.MaxSigmaValue = 3.5;

            var results = SpectralMerging.CombineSpectra(data.xArrays, data.yArrays, 
                tics, 60, options);

            double[] yAveraged = data.yArrays.CalculateAverage();

            

            double ticOriginal = data.yArrays[29].Sum();
            double[] originalDataNorm = data.yArrays[29].Select(i => i / ticOriginal).ToArray();

            MRSNoiseEstimator.MRSNoiseEstimation(yAveraged, 0.01, out double simpleAverage);
            MRSNoiseEstimator.MRSNoiseEstimation(results[1], 0.01, out double integration);
            MRSNoiseEstimator.MRSNoiseEstimation(originalDataNorm, 0.01, out double originalNoise);

            double integrationScale = Math.Sqrt(BiweightMidvariance(results[1]));
            double simpleAverageScale = Math.Sqrt(BiweightMidvariance(yAveraged));
            double originalScale = Math.Sqrt(BiweightMidvariance(originalDataNorm)); 

            double enr = (originalScale * integration) / (integrationScale * originalNoise);
            double simpleEnr = (originalScale * simpleAverage) / (simpleAverageScale * originalNoise);
            //peak.Plot().Show();
            //originalDataNorm.Plot(data.xArrays[0]).Show();
            //results.Plot().Show();
            //yAveraged.Plot(data.xArrays[0]).Show();

            Console.WriteLine("Original data to integrated ENR: {0}", enr);
            Console.WriteLine("Original signal to simple averaged ENR: {0}", simpleEnr);
            
            Console.WriteLine("Original Noise: {0} \nIntegration Noise: {1}", 
                originalNoise, integration);

            // get the noise mean 
            double snr = results[1][390..410].Average() / results[1][..150].StandardDeviation(); 
            double snrOriginal = originalDataNorm[390..410].Average() / originalDataNorm[..150].StandardDeviation();
            double snrImprovementSimple = yAveraged[390..410].Average() / yAveraged[..150].StandardDeviation();
            double improvement = (snr / snrOriginal);
            Console.WriteLine("Original SNR: {0}. Integrated SNR: {1}. \n SNR Improvement: {2}",
                snrOriginal, snr, improvement); 

            Console.WriteLine("Original SNR: {0}. Simple Average SNR: {1}. \n SNR Improvement: {2}", 
                snrOriginal, snrImprovementSimple, (snrImprovementSimple / snrOriginal));
            Console.WriteLine();
        }

        [Test, Repeat(1)]
        public void TestOptimizeSigmaMinMax()
        {
            ChromatographicPeak peak = new(30, 5, 60,
                0, null);
            Normal shotNoiseDistr = new(1, 0.05);
            Normal thermal = new(8000, 1200);
            TwoPeakTestData testData = new(60, 1000,
                (400, 10), (600, 10), peak, shotNoiseDistr);
            testData.CreateDataSet(thermal, 100000, 200000);

            var data = testData.GetJaggedArrays();
            var tics = testData.CalculateTics();
            double sigmaMin = 1.0;
            double sigmaMax = 1.0;

            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.PerformNormalization = true;
            options.BinSize = 1.0;
            options.RejectionType = RejectionType.SigmaClipping;
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.MaxSigmaValue = 1.5;
            options.MinSigmaValue = 1.5;

            double[] yAveraged = data.yArrays.CalculateAverage();

            double ticOriginal = data.yArrays[29].Sum();
            
            double[] originalDataNorm = data.yArrays[29].Select(i => i / ticOriginal).ToArray();
            double originalScale = Math.Sqrt(BiweightMidvariance(originalDataNorm));
            MRSNoiseEstimator.MRSNoiseEstimation(originalDataNorm, 0.01, out double originalNoise);

            double enr = 0;
            double nextEnr = 0;
            int interations = 0;
            double testValue = 0;
            double[][] results = new double[2][];
            do
            {
                options.MaxSigmaValue += 0.1;
                options.MinSigmaValue += 0.1;
                enr = nextEnr; 
                results = SpectralMerging.CombineSpectra(data.xArrays, data.yArrays,
                    tics, 60, options);
                MRSNoiseEstimator.MRSNoiseEstimation(results[1], 0.01, out double integration);
                double integrationScale = Math.Sqrt(BiweightMidvariance(results[1]));

                nextEnr = (originalScale * integration) / (integrationScale * originalNoise);

                testValue = (nextEnr - enr);
                Console.WriteLine("Iteration: {0}. Enr: {1}", interations, nextEnr); 
                interations++;

            } while (interations < 1000 && testValue > 0 && Math.Abs(testValue) > 0.0001);

            Console.WriteLine(options.MaxSigmaValue);
            Console.WriteLine(options.MinSigmaValue);
            Console.WriteLine(enr);
            results.Plot().Show();
        }

        [Test, TestCase(100)]
        public void TestMovingNoiseData(int numberSpectra)
        {
            // gaussian peak with new noise distribution that shifts x values differently each time 
            List<GaussianPeak> peakList = new();

            for (int i = 0; i < numberSpectra; i++)
            {
                Random random = new Random();
                GaussianPeak baseGaussian = new GaussianPeak(500, 50, 1000, 0, 1.0, null, 1000);
                // create a 50% chance that there will be a randomly occurring peak in the spectra. 
                if (random.NextDouble() < 0.25)
                {
                    int peakLocation = random.Next(250, 750);
                    int peakStddev = random.Next(1, 10);
                    int baseIntensity = random.Next(250, 750);

                    GaussianPeak noisePeak = new((double)peakLocation, (double)peakStddev, 1000, 0, 1.0, null,
                        intensity: (double)baseIntensity);

                    baseGaussian = (GaussianPeak)(baseGaussian + noisePeak); 
                }

                Normal xArrayNoise = new Normal(0, 0.01); 
                NoiseData noise = new NoiseData(1000, 0, 1.0d, 10, 1); 
                //noise.AddNoiseToXArray(xArrayNoise);
                peakList.Add((GaussianPeak)(baseGaussian + noise));
            }

            var data = peakList.GetJaggedArrays();

            var tics = peakList.CalculateTics();
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.PerformNormalization = true;
            options.BinSize = 1.0;
            options.RejectionType = RejectionType.WinsorizedSigmaClipping;
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.MinSigmaValue = 3.5;
            options.MaxSigmaValue = 3.5;


            var results = SpectralMerging.CombineSpectra(data.xArrays, data.yArrays,
                tics, numSpectra:numberSpectra, options);
            var averagedArray = data.yArrays.AverageJaggedArray();


            data.yArrays[0].Plot(data.xArrays[0]).Show();
            
            var plot = results.Plot();
            Layout layout = new Layout(); 
            layout.SetValue("title", "Winsorized sigma clipping (100 spectra).");
            layout.SetValue("plot_bgcolor", "#FFFFFF");
            plot.WithLayout(layout); 
            plot.Show();

            averagedArray.Plot(data.xArrays[0]).Show(); 

        }

        [Test, TestCase(100)]
        public void TestVaryingNoise(int numberSpectra)
        {
            // gaussian peak with new noise distribution that shifts x values differently each time 
            List<GaussianPeak> peakList = new();

            for (int i = 0; i < numberSpectra; i++)
            {
                Random random = new Random();
                GaussianPeak baseGaussian = new GaussianPeak(500, 50, 1000, 
                    0, 1.0, null, 100000);

                double noiseIntensity = (double)random.Next(250, 750);
                double noiseStddev = 0.25 * noiseIntensity;
                NoiseData noise = new NoiseData(1000, 0, 1.0d,
                    noiseIntensity, noiseStddev, random); 
                
                baseGaussian = (GaussianPeak)(baseGaussian + noise);
                
                //noise.AddNoiseToXArray(xArrayNoise);
                peakList.Add((GaussianPeak)(baseGaussian + noise));
            }

            var data = peakList.GetJaggedArrays();

            var tics = peakList.CalculateTics();
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.PerformNormalization = true;
            options.BinSize = 1.0;
            options.RejectionType = RejectionType.WinsorizedSigmaClipping;
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.MinSigmaValue = 3.5;
            options.MaxSigmaValue = 3.5;


            var results = SpectralMerging.CombineSpectra(data.xArrays, data.yArrays,
                tics, numSpectra: numberSpectra, options);
            var averagedArray = data.yArrays.AverageJaggedArray();


            data.yArrays[0].Plot(data.xArrays[0]).Show();

            var plot = results.Plot();
            Layout layout = new Layout();
            layout.SetValue("title", "Winsorized sigma clipping (100 spectra).");
            layout.SetValue("plot_bgcolor", "#FFFFFF");
            plot.WithLayout(layout);
            plot.Show();

            averagedArray.Plot(data.xArrays[0]).Show();

        }

        [Test]
        [TestCase(100)]
        public void TestAdditionOfManyShotPeaks(int numberSpectra)
        {
            // gaussian peak with new noise distribution that shifts x values differently each time 
            List<GaussianPeak> peakList = new();


            for (int i = 0; i < numberSpectra; i++)
            {
                Random random = new Random();
                GaussianPeak baseGaussian = new GaussianPeak(500, 50, 1000, 0, 1.0, null, 10000);
                
                // generate a random number of peaks 
                int numberOfRandomPeaks = random.Next(10, 50);
                for (int j = 0; j < numberOfRandomPeaks; j++)
                {
                    Random randomInnerLoop = new(); 
                    double peakLocation = randomInnerLoop.Next(100, 900);
                    double peakWidth = randomInnerLoop.Next(1, 5);
                    double peakIntensity = randomInnerLoop.Next(100, 400);
                    GaussianPeak peak = new(peakLocation, peakWidth, 1000, 0, 1.0, null, peakIntensity);
                    baseGaussian = (GaussianPeak)(baseGaussian + peak); 
                }

                NoiseData noise = new NoiseData(1000, 0, 1.0d, 10, 1);
                //noise.AddNoiseToXArray(xArrayNoise);
                peakList.Add((GaussianPeak)(baseGaussian + noise));
            }

            var data = peakList.GetJaggedArrays();

            var tics = peakList.CalculateTics();
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.PerformNormalization = true;
            options.BinSize = 1.0;
            options.RejectionType = RejectionType.WinsorizedSigmaClipping;
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.MinSigmaValue = 3.5;
            options.MaxSigmaValue = 3.5;


            var results = SpectralMerging.CombineSpectra(data.xArrays, data.yArrays,
                tics, numSpectra: numberSpectra, options);
            var averagedArray = data.yArrays.AverageJaggedArray();

            data.yArrays[0].Plot(data.xArrays[0]).Show();
            results.Plot().Show();
            averagedArray.Plot(data.xArrays[0]).Show();
        }



        private double MedianAbsoluteDeviationFromMedian(double[] array)
        {
            double arrayMedian = BasicStatistics.CalculateMedian(array);
            double[] results = new double[array.Length];
            for (int j = 0; j < array.Length; j++)
            {
                results[j] = Math.Abs(array[j] - arrayMedian);
            }

            return BasicStatistics.CalculateMedian(results);
        }
        private double BiweightMidvariance(double[] array)
        {
            double[] y_i = new double[array.Length];
            double[] a_i = new double[array.Length];
            double MAD_X = MedianAbsoluteDeviationFromMedian(array);
            double median = BasicStatistics.CalculateMedian(array);
            for (int i = 0; i < y_i.Length; i++)
            {
                y_i[i] = (array[i] - median) / (9d * MAD_X);
                if (y_i[i] < 1d)
                {
                    a_i[i] = 1d;
                }
                else
                {
                    a_i[i] = 0;
                }
            }

            // biweight midvariance calculation

            double denomSum = 0;
            double numeratorSum = 0;
            for (int i = 0; i < y_i.Length; i++)
            {
                numeratorSum += a_i[i] * Math.Pow(array[i] - median, 2) * Math.Pow(1 - y_i[i] * y_i[i], 4);
                denomSum += a_i[i] * (1 - 5 * y_i[i] * y_i[i]) * (1 - y_i[i] * y_i[i]);
            }

            return (double)y_i.Length * numeratorSum / Math.Pow(Math.Abs(denomSum), 2);
        }
    }

    


    public static class JaggedPlotter
    {
        public static double[] AverageJaggedArray(this double[][] yArrays)
        {
            double[] results = new double[yArrays[0].Length]; 
            for (int i = 0; i < yArrays.Length; i++)
            {
                double tic = yArrays[i].Sum(); 
                for (int j = 0; j < yArrays[i].Length; j++)
                {
                    results[j] += yArrays[i][j] / tic; 
                }
            }

            return results.Select(i => i / yArrays.Length).ToArray(); 
        }
        public static double[] CalculateTics(this List<GaussianPeak> peakList)
        {
            return peakList.Select(i => i.Yarray.Sum()).ToArray();
        }
        public static (double[][] xArrays, double[][] yArrays) GetJaggedArrays(this List<GaussianPeak> peakList)
        {
            double[][] xArrays = new double[peakList.Count][];
            double[][] yArrays = new double[peakList.Count][];

            for (int i = 0; i < peakList.Count; i++)
            {
                xArrays[i] = peakList[i].Xarray;
                yArrays[i] = peakList[i].Yarray;
            }
            return (xArrays, yArrays);
        }
        public static GenericChart.GenericChart Plot(this double[][] resultsJagged)
        {
            return Chart2D.Chart.Scatter<double, double, int>(resultsJagged[0], 
                resultsJagged[1], StyleParam.Mode.Lines);
        }

        public static double[] CalculateAverage(this double[][] yArrays)
        {
            double[] results = new double[yArrays[0].Length]; 
            for (int i = 0; i < yArrays.Length; i++)
            {
                double tic = yArrays[i].Sum(); 
                for (int j = 0; j < yArrays[i].Length; j++)
                {
                    results[j] += (yArrays[i][j] / tic); 
                }
            }

            return results.Select(i => i / yArrays.Length).ToArray(); 
        }

        public static GenericChart.GenericChart Plot(this double[] yArray, double[] xArray)
        {
            return Chart2D.Chart.Scatter<double, double, int>(xArray,
                yArray, StyleParam.Mode.Lines);
        }
        public static IList<double> FindPeaks(IList<double> values, int rangeOfPeaks)
        {
            List<double> peaks = new List<double>();

            int checksOnEachSide = rangeOfPeaks / 2;
            for (int i = 0; i < values.Count; i++)
            {
                double current = values[i];
                IEnumerable<double> range = values;
                if (i > checksOnEachSide)
                    range = range.Skip(i - checksOnEachSide);
                range = range.Take(rangeOfPeaks);
                if (current == range.Max())
                    peaks.Add(current);
            }
            return peaks;
        }

        public static double FindPeak1ToPeak2Ratio(this double[] yArray, int index1, int index2)
        {
            //var peaks = FindPeaks(yArray.ToList(), 100).ToList();
            //peaks.Sort();
            //var largestPeaks = peaks.TakeLast(2).ToArray(); 

            //return largestPeaks[0] / largestPeaks[1]; 
            return yArray[index1] / yArray[index2];
        }
    }
    public abstract class TestData
    {
        public double[] Xarray { get; set; }
        public double[] Yarray { get; set; }
        private int Length { get; set; }
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
        public (double,double) Peak1 { get; }
        public (double,double) Peak2 { get; }
        public List<(double Peak1, double Peak2, double ChromRatio)> PeakRatios { get; set; }
        public int Length { get; }
        public TwoPeakTestData(int numberSpectra, int length, (double, double) peak1, 
            (double,double) peak2, ChromatographicPeak chromatogram, 
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

        public void AddNoise(Normal noiseDistribution)
        {
            NoiseDistribution = noiseDistribution;
            AddNoise();
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
}
