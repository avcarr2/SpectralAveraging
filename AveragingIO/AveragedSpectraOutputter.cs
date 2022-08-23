using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using AveragingIO;
using IO.MzML;
using MassSpectrometry;
using Nett;

namespace SpectralAveraging
{
    public static class AveragedSpectraOutputter
    {
        public static void OutputAveragedScans(MsDataScan[] averagedScans, SpectralAveragingOptions options,
           string spectraPath)
        {
            switch (options.OutputType)
            {
                case OutputType.mzML:
                    OutputAveragedSpectraAsMzML(averagedScans, options, spectraPath);
                    break;

                case OutputType.txt:
                    OutputAveragedSpectraAsTxtFile(averagedScans, options, spectraPath);
                    break;


                default: throw new NotImplementedException("Output type not implemented");
            }

            if (options.OutputOptions)
            {
                string directoryPath = Path.GetDirectoryName(spectraPath);
                string optionsPath = Path.Combine(directoryPath, options.ToString() + ".toml");
                Toml.WriteFile(options, optionsPath);
            }
        }

        private static void OutputAveragedSpectraAsMzML(MsDataScan[] averagedScans, SpectralAveragingOptions options,
            string spectraPath)
        {
            SourceFile sourceFile = SpectraFileHandler.GetSourceFile(spectraPath);
            MsDataFile msDataFile = new(averagedScans, sourceFile);
            string spectraDirectory = Path.GetDirectoryName(spectraPath);
            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + Path.GetFileNameWithoutExtension(spectraPath) + ".mzML");

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, averagedPath, true);
        }

        private static void OutputAveragedSpectraAsTxtFile(MsDataScan[] averagedScans, SpectralAveragingOptions options,
            string spectraPath)
        {
            string spectraDirectory = Path.GetDirectoryName(spectraPath);
            if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
            {
                spectraDirectory = Path.Combine(spectraDirectory, "AveragedSpectra");
                if (!Directory.Exists(spectraDirectory))
                    Directory.CreateDirectory(spectraDirectory);
            }

            string averagedPath = Path.Combine(spectraDirectory,
                "Averaged_" + Path.GetFileNameWithoutExtension(spectraPath) + options + ".txt");

            foreach (var scan in averagedScans)
            {
                if (options.SpectraFileProcessingType != SpectraFileProcessingType.AverageAll)
                    averagedPath = Path.Combine(spectraDirectory,
                        "Averaged_" + Path.GetFileNameWithoutExtension(spectraPath) + "_" + scan.OneBasedScanNumber + options + ".txt");
                using (StreamWriter writer = new StreamWriter(File.Create(averagedPath)))
                {

                    for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                    {
                        writer.WriteLine(scan.MassSpectrum.XArray[i] + "\t" + scan.MassSpectrum.YArray[i]);
                    }
                }
            }
        }
    }
}
