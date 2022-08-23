using System.Text;

namespace SpectralAveraging
{
    public class SpectralAveragingOptions 
    {
        #region Averaging Options
        public RejectionType RejectionType { get; set; }
        public WeightingType WeightingType { get; set; }
        public SpectrumMergingType SpectrumMergingType { get; set; }
        public bool PerformNormalization { get; set; }
        public double Percentile { get; set; }
        public double MinSigmaValue { get; set; }
        public double MaxSigmaValue { get; set; }
        public double BinSize { get; set; }

        #endregion

        #region File Processing Options
        public SpectraFileProcessingType SpectraFileProcessingType { get; set; }
        public int NumberOfScansToAverage { get; set; }
        public int ScanOverlap { get; set; }
        public OutputType OutputType { get; set; }
        public bool OutputOptions { get; set; }

        #endregion


        public SpectralAveragingOptions()
        {
            SetDefaultValues();
        }

        /// <summary>
        /// Can be used to set the values of the options class in one method call
        /// </summary>
        /// <param name="rejectionType">rejection type to be used</param>
        /// <param name="percentile">percentile for percentile clipping rejection type</param>
        /// <param name="sigma">sigma value for sigma clipping rejection types</param>
        public void SetValues(RejectionType rejectionType = RejectionType.NoRejection,
            WeightingType intensityWeighingType = WeightingType.NoWeight, SpectrumMergingType spectrumMergingType = SpectrumMergingType.SpectrumBinning,
            bool performNormalization = true, double percentile = 0.1, double minSigma = 1.5, double maxSigma = 1.5, double binSize = 0.01,
            SpectraFileProcessingType spectraFileProcessingType = SpectraFileProcessingType.AverageAll, int numberOfScansToAverage = 5, 
            int scanOverlap = 2, OutputType outputType = OutputType.mzML, bool outputOptions = false)
        {
            RejectionType = rejectionType;
            WeightingType = intensityWeighingType;
            SpectrumMergingType = spectrumMergingType;
            PerformNormalization = performNormalization;
            Percentile = percentile;
            MinSigmaValue = minSigma;
            MaxSigmaValue = maxSigma;
            BinSize = binSize;
            SpectraFileProcessingType = spectraFileProcessingType;
            NumberOfScansToAverage = numberOfScansToAverage;
            ScanOverlap = scanOverlap;
            OutputType = outputType;
            OutputOptions = outputOptions;
        }

        /// <summary>
        /// Sets the values of the options to their defaults
        /// </summary>
        public void SetDefaultValues()
        {
            RejectionType = RejectionType.NoRejection;
            WeightingType = WeightingType.NoWeight;
            SpectrumMergingType = SpectrumMergingType.SpectrumBinning;
            PerformNormalization = true;
            Percentile = 0.1;
            MinSigmaValue = 1.5;
            MaxSigmaValue = 1.5;
            BinSize = 0.01;
            SpectraFileProcessingType = SpectraFileProcessingType.AverageAll;
            NumberOfScansToAverage = 5;
            ScanOverlap = 2;
            OutputType = OutputType.mzML;
            OutputOptions = false;
        }

        /// <summary>
        /// Override for the ToString method that can be used for file output naming
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.Append(RejectionType.ToString() + '_');
            stringBuilder.Append(WeightingType.ToString() + '_');
            if (PerformNormalization)
                stringBuilder.Append("Normalized_");

            // rejection type specific 
            if (RejectionType == RejectionType.PercentileClipping)
                stringBuilder.Append("Percentile-" + Percentile + '_');
            if (RejectionType is RejectionType.WinsorizedSigmaClipping or RejectionType.AveragedSigmaClipping
                or RejectionType.SigmaClipping)
            {
                stringBuilder.Append("MinSigma-" + MinSigmaValue + '_');
                stringBuilder.Append("MaxSigma-" + MaxSigmaValue + '_');
            }

            stringBuilder.Append("BinSize-" + BinSize + '_');

            // file processing specific
            stringBuilder.Append(SpectraFileProcessingType.ToString() + '_');
            if (SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
            {
                stringBuilder.Append("Averaged" + NumberOfScansToAverage + "Scans_");
                if (SpectraFileProcessingType.ToString().Contains("Overlap"))
                {
                    stringBuilder.Append("ScanOverlap-" + ScanOverlap);
                }
            }

            return stringBuilder.ToString();
        }
    }


}
