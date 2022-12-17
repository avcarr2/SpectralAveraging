using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Nett;
using SpectralAveraging.NoiseEstimates;

namespace SpectralAveraging.DataStructures
{

    public class BinnedSpectra
    {
        public List<PixelStack> PixelStacks { get; set; }
        public Dictionary<int, double> NoiseEstimates { get; private set; }
        public Dictionary<int, double> ScaleEstimates { get; private set; }
        public Dictionary<int, double> Weights { get; private set; }
        private int ReferenceSpectraIndex { get; set; }
        public double[] Tics { get; private set; }
        public int NumSpectra => PixelStacks[0].Intensity.Count;
        private List<double[]> RecalculatedSpectra { get; set; }

        public BinnedSpectra(int referenceSpectraIndex = 0)
        {
            PixelStacks = new List<PixelStack>();
            NoiseEstimates = new Dictionary<int, double>();
            ScaleEstimates = new Dictionary<int, double>();
            Weights = new Dictionary<int, double>();
            ReferenceSpectraIndex = referenceSpectraIndex; 
        }

        public void ProcessPixelStacks(SpectralAveragingOptions options)
        {
            Parallel.ForEach(PixelStacks, pixelStack =>
            {
                pixelStack.PerformRejection(options);
            }); 
        }

        public void ConsumeSpectra(double[][] xArrays, double[][] yArrays,
            int numSpectra, double binSize)
        {
            double min = 100000;
            double max = 0;
            for (int i = 0; i < numSpectra; i++)
            {
                min = Math.Min(xArrays[i][0], min);
                max = Math.Max(xArrays[i].Max(), max);
            }

            int numberOfBins = (int)Math.Ceiling((max - min) * (1 / binSize));
            // go through each scan and place each (m/z, int) from the spectra into a jagged array

            // 1) find all values of x that fall within a bin.

            List<List<BinValue>> listBinValues = new();

            for (int i = 0; i < numSpectra; i++)
            {
                // currently this updates the y value with the most recent value from the array. 
                // what it really needs to do is a linear interpolation. 

                listBinValues.Add(CreateBinValueList(xArrays[i], yArrays[i], min, binSize));
            }

            int[] spectraId = Enumerable.Range(0, numSpectra).ToArray(); 

            for (int i = 0; i < numberOfBins; i++)
            {
                List<double> xVals = new();
                List<double> yVals = new();
                foreach (List<BinValue> binValList in listBinValues)
                {
                    var valsInBin = binValList
                        .Where(m => m.Bin == i);
                    if (valsInBin.Count() > 1)
                    {
                        xVals.Add(valsInBin.Average(m => m.Mz));
                        yVals.Add(valsInBin.Average(m => m.Intensity));
                    }
                    else
                    {
                        xVals.Add(valsInBin.First().Mz);
                        yVals.Add(valsInBin.First().Intensity);
                    }
                }
                PixelStacks.Add(new PixelStack(xVals, yVals, spectraId));
            }
        }

        public void PerformNormalization()
        {
            Parallel.For(0, PixelStacks.Count, i =>
            {
                // pixelStacks[i].Length will be the number of spectra to be averaged, 
                // which is equal to tics. 
                for (int j = 0; j < PixelStacks[i].Length; j++)
                {
                    PixelStacks[i].Intensity[j] /= Tics[j];
                }
            }); 
        }

        public void CalculateNoiseEstimates(WaveletType waveletType = WaveletType.Haar, 
            double epsilon = 0.01, int maxIterations = 25)
        {
            // note that noise estimate needs to be run on the pixel stacks. 
            for (int i = 0; i < NumSpectra; i++)
            {
                double[] tempValArray = PopIntensityValuesFromPixelStackList(i);
                bool success = MRSNoiseEstimator.MRSNoiseEstimation(tempValArray, epsilon, out double noiseEstimate,
                    waveletType: waveletType, maxIterations: maxIterations);
                if (!success || double.IsNaN(noiseEstimate))
                {
                    noiseEstimate = BasicStatistics.CalculateStandardDeviation(tempValArray); 
                }
                NoiseEstimates.Add(i, noiseEstimate);
            }
        }

        public void CalculateScaleEstimates()
        {
            double reference = 0;
            for (int i = 0; i < NumSpectra; i++)
            {
                double[] tempValArray = PopIntensityValuesFromPixelStackList(i);
                double scale = BiweightMidvariance(tempValArray);
                if (i == 0)
                {
                    reference = scale;
                }
                ScaleEstimates.Add(i, reference / scale);
            }
            
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

        public void CalculateWeights()
        {
            double referenceScale = Weights[ReferenceSpectraIndex]; 
            foreach (var entry in NoiseEstimates)
            {
                var successScale = ScaleEstimates.TryGetValue(entry.Key,
                    out double scale);
                if (!successScale) continue;

                var successNoise = NoiseEstimates.TryGetValue(entry.Key,
                    out double noise);
                if (!successNoise) continue;

                double scaleFactor = referenceScale / scale;
                double weight = 1d / Math.Pow((scaleFactor * noise), 2);

                Weights.Add(entry.Key, weight);
            }
        }

        public void RecalculateTics()
        {
            Tics = new double[NumSpectra];
            RecalculatedSpectra = new List<double[]>(); 
            for (int i = 0; i < NumSpectra; i++)
            {
                RecalculatedSpectra.Add(PopIntensityValuesFromPixelStackList(i));
                Tics[i] = RecalculatedSpectra[i].Sum();
            }
        }

        public void MergeSpectra()
        {
            var weights = Weights.Values.ToArray(); 
            Parallel.ForEach(PixelStacks, pixelStack =>
            {
                pixelStack.Average(weights);
            });
        }

        public double[][] GetMergedSpectrum()
        {
            double[] xArray = PixelStacks.Select(i => i.Mz).ToArray();
            double[] yArray = PixelStacks.Select(i => i.MergedValue).ToArray();
            return new[] { xArray, yArray };
        }
        private static List<BinValue> CreateBinValueList(double[] xArray, double[] yArray,
            double min, double binSize)
        {
            var binIndices = xArray
                .Select(i => (int)Math.Floor((i - min) / binSize))
                .ToArray();
            List<BinValue> binValues = new List<BinValue>();

            for (int i = 0; i < binIndices.Length; i++)
            {
                binValues.Add(new BinValue(binIndices[i], xArray[i], yArray[i]));
            }
            return binValues;
        }

        private double[] PopIntensityValuesFromPixelStackList(int index)
        {
            double[] results = new double[PixelStacks.Count];
            for (int i = 0; i < PixelStacks.Count; i++)
            {
                results[i] = PixelStacks[i].Intensity[index]; 
            }

            return results; 
        }
    }
    /// <summary>
    /// Record type used to facilitate bin, mz, and intensity matching. 
    /// </summary>
    /// <param name="Bin"></param>
    /// <param name="Mz"></param>
    /// <param name="Intensity"></param>
    internal readonly record struct BinValue(int Bin, double Mz, double Intensity);
}
