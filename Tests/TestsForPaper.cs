using System.Xml.Schema;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.Statistics;
using Microsoft.VisualStudio.TestPlatform.CommunicationUtilities;
using NUnit.Framework.Interfaces;
using Plotly.NET;
using Plotly.NET.LayoutObjects;
using SpectralAveraging;
using SpectralAveraging.NoiseEstimates;
using MzLibUtil;
using Chemistry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using MassSpectrometry;

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

        //[Test]
        //[TestCase(100)]
        //public void TestAdditionOfManyShotPeaks(int numberSpectra)
        //{
        //    // gaussian peak with new noise distribution that shifts x values differently each time 
        //    var peakList = CreateListGaussianPeaksWithNoise(numberSpectra, 500, 50, 1000, 1.0);

        //    var data = peakList.GetJaggedArrays();

        //    var tics = peakList.CalculateTics();
        //    var winOptions = CreateWinsorizedSigmaClippingOptions(); 

        //    var results = SpectralMerging.CombineSpectra(data.xArrays, data.yArrays,
        //        tics, numSpectra: numberSpectra, winOptions);
        //    var averagedArray = data.yArrays.AverageJaggedArray();

        //    data.yArrays[0].Plot(data.xArrays[0]).Show();
        //    results.Plot().Show();
        //    averagedArray.Plot(data.xArrays[0]).Show();
        //}

        private static void AddLayoutToChart()
        {

        }

        private static double CalculateENR(double[] referenceArray, double[] testArray,
            int lowIndex, int highIndex)
        {
            double medianTestArray = BasicStatistics.CalculateMedian(testArray[lowIndex..highIndex]); 
            double medianRefArray = BasicStatistics.CalculateMedian(referenceArray[lowIndex..highIndex]);

            double diff = medianTestArray - medianRefArray;

            referenceArray = referenceArray.Select(i => i - diff).ToArray(); 

            double referenceScaleEstimate = EstimateScale(referenceArray, lowIndex, highIndex);
            double testScaleEstimate = EstimateScale(testArray, lowIndex, highIndex); 
            
            double k = referenceScaleEstimate / testScaleEstimate;

            MRSNoiseEstimator.MRSNoiseEstimation(referenceArray[lowIndex..highIndex], 0.01, out double refNoiseEstimate);
            MRSNoiseEstimator.MRSNoiseEstimation(testArray[lowIndex..highIndex], 0.01, out double testNoiseEstimate);

            return (k * testNoiseEstimate) / refNoiseEstimate;
        }

        private static double EstimateScale(double[] array, int lowIndex, int highIndex)
        {
            return Math.Sqrt(BiweightMidvariance(array[lowIndex..highIndex]));
        }

        private static SpectralAveragingOptions CreateWinsorizedSigmaClippingOptions(double binSize)
        {
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.PerformNormalization = true;
            options.BinSize = binSize;
            options.RejectionType = RejectionType.WinsorizedSigmaClipping;
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.MinSigmaValue = 3.5;
            options.MaxSigmaValue = 3.5;
            
            return options;
        }

        //private static List<GaussianPeak> CreateListGaussianPeaksWithNoise(int numberSpectra,
        //    double meanBase, double stddevBase, double intensityBase,
        //    double peakLowVal = 100, double peakHighVal = 900, double peakWidthLow = 1,
        //    double peakWidthHigh = 5, int peakIntLow = 100, int peakIntHigh = 400,
        //    int length = 1000, int startVal = 0, double spacing = 1.0)
        //{
        //    List<GaussianPeak> peakList = new();
        //    for (int i = 0; i < numberSpectra; i++)
        //    {
        //        peakList.Add(GenerateGuassianPeakWithRandomNoise(meanBase, stddevBase, 
        //            intensityBase, peakLowVal, peakHighVal, peakWidthLow, peakWidthHigh,
        //            peakIntLow, peakIntHigh, length, startVal, spacing));
        //    }
            
        //    return peakList;
        //}

        
        private static double MedianAbsoluteDeviationFromMedian(double[] array)
        {
            double arrayMedian = BasicStatistics.CalculateMedian(array);
            double[] results = new double[array.Length];
            for (int j = 0; j < array.Length; j++)
            {
                results[j] = Math.Abs(array[j] - arrayMedian);
            }

            return BasicStatistics.CalculateMedian(results);
        }
        private static double BiweightMidvariance(double[] array)
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
        

        [Test]
        [TestCase(25)]
        [Repeat(1)]
        public void Figure1(int numberSpectra)
        {
            // Need: 
            /*
             *  1) SNR Estimates
             *  2) ENR Estimates
             *  3) Three charts in a row: Representative, Averaging with no rejection, Averaging with rejection.
             *
             */

            // generate a list of spectra. 
            List<GaussianPeak> peakList = new();
            double[] tics = new double[numberSpectra];
            for (int i = 0; i < numberSpectra; i++)
            {
                // create a spectra and add shot peaks to it: 
                // first create the noiseless distribution: 
                GaussianPeak peak = new GaussianPeak(500, 15, 1000, 
                    0.0, 1.0, null, 100000); 
                
                // add high frequency noise distribution 
                peak.AddNoise(new Normal(1000, 50));

                // add a random amount of low-frequency peaks 
                peak.AddShotPeaks(5, 50, 200, 
                    800, 1, 3, 100, 20000);
                tics[i] = peak.Yarray.Sum(); 
                peakList.Add(peak);
            }

            // perform naive averaging 
            (double[][] xArrays, double[][] yArrays) arrays = peakList.GetJaggedArrays();
            double[] naiveAverage = arrays.yArrays.CalculateAverage();
            // perform averaging with rejection


            // calculate the ENR estimate for the 200-800 m/z range 
            var options = CreateWinsorizedSigmaClippingOptions(0.5);
            double[][] outputArrays = SpectralMerging.CombineSpectra(arrays.xArrays, arrays.yArrays, tics, 25,
                options);

            // calculate the ENR estiamte for the 0-200 m/z range 
            // ENR estimates for the averaging with rejection
            double[] normalizedRefArray = arrays.yArrays[0].Select(i => i / tics[0]).ToArray();

            double enrLHighreqNoise = CalculateENR(normalizedRefArray, outputArrays[1], 0, 199);
            double enrLowFreqNoise = CalculateENR(normalizedRefArray, outputArrays[1], 200, 800);
            // ENR estimates for naive averaging 
            double enrHighFreqNoiseNaive = CalculateENR(normalizedRefArray, naiveAverage, 0, 199);
            double enrLowFreqNoiseNaive = CalculateENR(normalizedRefArray, naiveAverage, 200, 800);
            // ENR estimate between naive and with rejection 
            double enrHighFreqNoiseComp = CalculateENR(naiveAverage, outputArrays[1], 0, 199);
            double enrLowFreqNoiseComp = CalculateENR(naiveAverage, outputArrays[1], 200, 800);

        }

        [Test]
        [TestCase(25)]
        // simulated protein
        public void Figure2(int numberSpectra)
        {
            UsefulProteomicsDatabases.PeriodicTableLoader.Load();
            Peptide protein = new Peptide("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
            var distribution = IsotopicDistribution.GetDistribution(protein.GetChemicalFormula());
            var transformedIntensities = distribution.Intensities.Select(i => i * 1000);
            List<SimulatedProtein> peakList = new();
            double[] tics = new double[numberSpectra];
            double max = distribution.Intensities.Max(); 
            for (int i = 0; i < numberSpectra; i++)
            {
                // create a spectra and add shot peaks to it: 
                // first create the noiseless distribution: 
                SimulatedProtein peak = new SimulatedProtein(distribution, 2500, 0.01);

                // add high frequency noise distribution 
                peak.AddNoise(new Normal(0.1, 0.01));

                // add a random amount of low-frequency peaks 
                peak.AddShotPeaks(5, 50, 
                    peak.Xarray.Min(),
                    peak.Xarray.Max(),
                    0.1, 0.01, max/2, max);
                tics[i] = peak.Yarray.Sum();
                peakList.Add(peak);
            }
            (double[][] xArrays, double[][] yArrays) arrays = peakList.GetJaggedArrays();
            var options = CreateWinsorizedSigmaClippingOptions(0.01);
            double[][] outputArrays = SpectralMerging.CombineSpectra(arrays.xArrays, arrays.yArrays, 
                tics, 25, options);
            double[] naiveAverage = arrays.yArrays.CalculateAverage();
            double[] normalizedRefArray = arrays.yArrays[0].Select(i => i / tics[0]).ToArray();

            normalizedRefArray.Plot(arrays.xArrays[0]).Show();
            naiveAverage.Plot(arrays.xArrays[0]).Show(); 
            outputArrays.Plot().Show();
        }

        [Test]
        public void DemonstrateSimulatedProtein()
        {
            UsefulProteomicsDatabases.PeriodicTableLoader.Load();
            Peptide protein = new Peptide("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" +
                                          "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
            var distribution = IsotopicDistribution.GetDistribution(protein.GetChemicalFormula());

            SimulatedProtein simulated = new(distribution, 2500, 0.01);
            simulated.Plot().Show();
            simulated.AddShotPeaks(5, 25,
                simulated.Xarray.Min(),
                simulated.Xarray.Max(), 0.01, 0.01,
                0.001, 0.01, 2500, simulated.Xarray.Min(), 0.01);
            simulated.Plot().Show();

        }
    }

}
