﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using AveragingIO;
using IO.MzML;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace SpectralAveraging
{
    /// <summary>
    /// Manages the spectra files being averaged
    /// </summary>
    public static class SpectraFileProcessing
    {
        public static MsDataScan[] ProcessSpectra(List<MsDataScan> scans, SpectralAveragingOptions options)
        {
            switch (options.SpectraFileProcessingType)
            {
                case SpectraFileProcessingType.AverageAll:
                    return AverageAll(scans, options);

                case SpectraFileProcessingType.AverageEverynScans:
                    options.ScanOverlap = 0;
                    return AverageEverynScans(scans, options);

                case SpectraFileProcessingType.AverageEverynScansWithOverlap:
                    return AverageEverynScans(scans, options);

                case SpectraFileProcessingType.AverageDDAScans:
                    options.ScanOverlap = 0;
                    return AverageDDAScans(scans, options);

                case SpectraFileProcessingType.AverageDDAScansWithOverlap:
                    return AverageDDAScans(scans, options);

                default: throw new ArgumentOutOfRangeException(nameof(options));
            }
        }

        private static MsDataScan[] AverageAll(List<MsDataScan> scans, SpectralAveragingOptions options)
        {
            // average spectrum
            MsDataScan representativeScan = scans.First();
            MultiScanDataObject multiScanDataObject = new(SingleScanDataObject.ConvertMSDataScansInBulk(scans));
            MzSpectrum averagedSpectrum = SpectralMerging.CombineSpectra(multiScanDataObject, options);

            // create output
            MsDataScan averagedScan = new(averagedSpectrum, 1, representativeScan.OneBasedScanNumber,
                representativeScan.IsCentroid, representativeScan.Polarity, scans.Select(p => p.RetentionTime).Average(),
                averagedSpectrum.Range, null, representativeScan.MzAnalyzer, (double)multiScanDataObject.AverageIonCurrent,
                representativeScan.InjectionTime, null, representativeScan.NativeId);
            MsDataScan[] msDataScans = new MsDataScan[] { averagedScan };
            return msDataScans;
        }

        private static MsDataScan[] AverageEverynScans(List<MsDataScan> scans, SpectralAveragingOptions options)
        {
            List<MsDataScan> averagedScans = new();
            int scanNumberIndex = 1;
            for (int i = 0; i < scans.Count; i += options.NumberOfScansToAverage - options.ScanOverlap)
            {
                // get the scans to be averaged
                List<MsDataScan> scansToProcess = new();
                if (i <= options.ScanOverlap) // very start of the file
                {
                    scansToProcess = scans.GetRange(i, options.NumberOfScansToAverage);
                }
                else if (i + options.NumberOfScansToAverage > scans.Count) // very end of the file
                {
                    break;
                }
                else // anywhere in the middle of the file
                {
                    scansToProcess = scans.GetRange(i , options.NumberOfScansToAverage);
                }

                // average scans
                MsDataScan representativeScan = scansToProcess.First();
                MultiScanDataObject multiScanDataObject = new(SingleScanDataObject.ConvertMSDataScansInBulk(scansToProcess));
                MzSpectrum averagedSpectrum = SpectralMerging.CombineSpectra(multiScanDataObject, options);
                MsDataScan averagedScan = new(averagedSpectrum, scanNumberIndex, 1,
                    representativeScan.IsCentroid, representativeScan.Polarity, scansToProcess.Select(p => p.RetentionTime).Minimum(),
                    averagedSpectrum.Range, null, representativeScan.MzAnalyzer, (double)multiScanDataObject.AverageIonCurrent,
                    scansToProcess.Select(p => p.InjectionTime).Average(), null, representativeScan.NativeId);

                averagedScans.Add(averagedScan);
                scanNumberIndex++;
            }

            return averagedScans.ToArray();
        }

        private static MsDataScan[] AverageDDAScans(List<MsDataScan> scans, SpectralAveragingOptions options)
        {
            List<MsDataScan> averagedScans = new();
            List<MsDataScan> ms1Scans = scans.Where(p => p.MsnOrder == 1).ToList();
            List<MsDataScan> ms2Scans = scans.Where(p => p.MsnOrder == 2).ToList();
            List<MsDataScan> scansToProcess = new();

            int scanNumberIndex = 1;
            for (int i = 0; i < ms1Scans.Count; i += options.NumberOfScansToAverage - options.ScanOverlap)
            {
                // get the scans to be averaged
                scansToProcess.Clear();
                if (i <= options.ScanOverlap) // very start of the file
                {
                    scansToProcess = ms1Scans.GetRange(i, options.NumberOfScansToAverage);
                }
                else if (i + options.NumberOfScansToAverage > ms1Scans.Count) // very end of the file
                {
                    break;
                }
                else // anywhere in the middle of the file
                {
                    scansToProcess = ms1Scans.GetRange(i, options.NumberOfScansToAverage);
                }

                // average scans and add to averaged list
                MsDataScan representativeScan = scansToProcess.First();
                MultiScanDataObject multiScanDataObject = new(SingleScanDataObject.ConvertMSDataScansInBulk(scansToProcess));
                MzSpectrum averagedSpectrum = SpectralMerging.CombineSpectra(multiScanDataObject, options);
                MsDataScan averagedScan = new(averagedSpectrum, scanNumberIndex, 1,
                    representativeScan.IsCentroid, representativeScan.Polarity, representativeScan.RetentionTime,
                    averagedSpectrum.Range, null, representativeScan.MzAnalyzer, (double)multiScanDataObject.AverageIonCurrent,
                    scansToProcess.Select(p => p.InjectionTime).Average(), representativeScan.NoiseData, 
                    representativeScan.NativeId, representativeScan.SelectedIonMZ, representativeScan.SelectedIonChargeStateGuess, 
                    representativeScan.SelectedIonIntensity, representativeScan.IsolationMz, representativeScan.IsolationWidth, 
                    representativeScan.DissociationType, representativeScan.OneBasedPrecursorScanNumber, representativeScan.SelectedIonMonoisotopicGuessIntensity);
                averagedScans.Add(averagedScan);
                int precursorScanIndex = scanNumberIndex;
                scanNumberIndex++;

                // add the respective Ms2 scans
                IEnumerable<MsDataScan> ms2ScansFromAveragedScans = ms2Scans.Where(p =>
                    scansToProcess.Any(m => m.OneBasedScanNumber == p.OneBasedPrecursorScanNumber));
                foreach (var scan in ms2ScansFromAveragedScans)
                {
                    MsDataScan newScan = new(scan.MassSpectrum, scanNumberIndex, scan.MsnOrder, scan.IsCentroid,
                        scan.Polarity, scan.RetentionTime, scan.ScanWindowRange, scan.ScanFilter, scan.MzAnalyzer,
                        scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId, scan.SelectedIonMZ, 
                        scan.SelectedIonChargeStateGuess, scan.SelectedIonIntensity, scan.IsolationMz, scan.IsolationWidth, 
                        scan.DissociationType, precursorScanIndex, scan.SelectedIonMonoisotopicGuessMz);
                    averagedScans.Add(newScan);
                    scanNumberIndex++;
                }
            }

            var test = averagedScans.Select(p => (p.OneBasedScanNumber, p.OneBasedPrecursorScanNumber));
            return averagedScans.ToArray();
        }
    }
}
