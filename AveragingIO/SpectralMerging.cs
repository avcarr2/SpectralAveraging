using System.Collections.Concurrent;
using MassSpectrometry;
using MzLibUtil;
using System.Linq;

namespace SpectralAveraging
{
    public static class SpectralMerging
    {
        /// <summary>
        /// Calls the specific merging function based upon the current static field SpecrumMergingType
        /// </summary>
        /// <param name="scans"></param>
        public static MzSpectrum CombineSpectra(double[][] xArrays, double[][] yArrays, double[] totalIonCurrents, int numSpectra, SpectralAveragingOptions options)
        {
            MzSpectrum compositeSpectrum = null;
            switch (options.SpectrumMergingType)
            {
                case SpectrumMergingType.SpectrumBinning:
                    compositeSpectrum = SpectrumBinning(xArrays, yArrays, totalIonCurrents, options.BinSize, numSpectra, options);
                    break;


                case SpectrumMergingType.MostSimilarSpectrum:
                    compositeSpectrum = MostSimilarSpectrum();
                    break;
            }
            return compositeSpectrum;
        }

        public static MzSpectrum CombineSpectra(MultiScanDataObject multiScanDataObject, SpectralAveragingOptions options)
        {
            return CombineSpectra(multiScanDataObject.XArrays, multiScanDataObject.YArrays, multiScanDataObject.TotalIonCurrent,
                multiScanDataObject.ScansToProcess, options);
        }

        /// <summary>
        /// Merges spectra into a two dimensional array of (m/z, int) values based upon their bin 
        /// </summary>
        /// <param name="scans">scans to be combined</param>
        /// <returns>MSDataScan with merged values</returns>
        public static MzSpectrum SpectrumBinning(double[][] xArrays, double[][] yArrays, double[] totalIonCurrents, double binSize, int numSpectra,
            SpectralAveragingOptions options)
        {
            // normalize if selected
            if (options.PerformNormalization)
            {
                for (int i = 0; i < xArrays.Length; i++)
                {
                    SpectrumNormalization.NormalizeSpectrumToTic(yArrays[i], totalIonCurrents[i], totalIonCurrents.Average());
                }
            }

            // calculate the bins to be utilized
            double min = 100000;
            double max = 0;
            for (int i = 0; i < numSpectra; i++)
            {
                min = Math.Min(xArrays[i][0], min);
                max = Math.Max(xArrays[i].Max(), max);
            }
            int numberOfBins = (int)Math.Ceiling((max - min) * (1 / binSize));

            double[][] xValuesBin = new double[numberOfBins][];
            double[][] yValuesBin = new double[numberOfBins][];
            // go through each scan and place each (m/z, int) from the spectra into a jagged array
            for (int i = 0; i < numSpectra; i++)
            {
                for (int j = 0; j < xArrays[i].Length; j++)
                {
                    int binIndex = (int)Math.Floor((xArrays[i][j] - min) / binSize);
                    if (xValuesBin[binIndex] == null)
                    {
                        xValuesBin[binIndex] = new double[numSpectra];
                        yValuesBin[binIndex] = new double[numSpectra];
                    }
                    xValuesBin[binIndex][i] = xArrays[i][j];
                    yValuesBin[binIndex][i] = yArrays[i][j];
                }
            }

            xValuesBin = xValuesBin.Where(p => p != null).ToArray();
            yValuesBin = yValuesBin.Where(p => p != null).ToArray();

            // average the remaining arrays to create the composite spectrum
            // this will clipping and averaging for y values as indicated in the settings
            double[] xArray = new double[xValuesBin.Length];
            double[] yArray = new double[yValuesBin.Length];

            for (int i = 0; i < yValuesBin.Length; i++)
            {
                xArray[i] = xValuesBin[i].Where(p => p != 0).Average();
            }

            var rangePartitioner = Partitioner.Create(0, yValuesBin.GetLength(0));
            Parallel.ForEach(rangePartitioner, (range, loopState) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    yArray[i] = ProcessSingleMzArray(yValuesBin[i], options); 
                }
            }); 


            // Create new MsDataScan to return
            MzRange range = new(min, max);
            MzSpectrum mergedSpectra = new(xArray, yArray, true);

            return mergedSpectra;
        }

        public static MzSpectrum MostSimilarSpectrum()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Main Engine of this Binning method, processes a single array of intesnity values for a single mz and returns their average
        /// </summary>
        /// <param name="intInitialArray"></param>
        /// <returns></returns>
        public static double ProcessSingleMzArray(double[] intInitialArray, SpectralAveragingOptions options)
        {
            double average;
            double[] weights;
            double[] trimmedArray;

            if (intInitialArray.Where(p => p != 0).Count() <= 1)
                return 0;

            trimmedArray = OutlierRejection.RejectOutliers(intInitialArray, options);
            if (trimmedArray.Where(p => p != 0).Count() <= 1)
                return 0;
            weights = BinWeighting.CalculateWeights(trimmedArray, options.WeightingType);
            average = MergePeakValuesToAverage(trimmedArray, weights);
            return average;
        }

        /// <summary>
        /// Calculates the weighted average value for each m/z point passed to it
        /// </summary>
        /// <param name="intValues">array of mz values to evaluate </param>
        /// <param name="weights">relative weight assigned to each of the mz values</param>
        /// <returns></returns>
        public static double MergePeakValuesToAverage(double[] intValues, double[] weights)
        {
            double numerator = 0;
            for (int i = 0; i < intValues.Count(); i++)
            {
                numerator += intValues[i] * weights[i];
            }
            double average = numerator / weights.Sum();
            return average;
        }
    }
}
