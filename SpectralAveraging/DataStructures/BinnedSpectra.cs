﻿using System;
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

        public BinnedSpectra()
        {
            PixelStacks = new List<PixelStack>();
            NoiseEstimates = new Dictionary<int, double>();
            ScaleEstimates = new Dictionary<int, double>();
            Weights = new Dictionary<int, double>();
        }

        public void ProcessPixelStacks(SpectralAveragingOptions options)
        {
            foreach (var pixelStack in PixelStacks)
            {
                pixelStack.PerformRejection(options);
            }
        }

        public void ConsumeSpectra(double[][] xArrays, double[][] yArrays, 
            double[] totalIonCurrents, int numSpectra, double binSize)
        {
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

            // 1) find all values of x that fall within a bin.
            
            // 2) perform a cubic spline interpolation to find the value of the 
            // central point (aka the bin center). 

            for (int i = 0; i < numSpectra; i++)
            {
                // currently this updates the y value with the most recent value from the array. 
                // what it really needs to do is a linear interpolation. 

                var binIndices = xArrays[i]
                    .Select(i => (int)Math.Floor((i - min) / binSize))
                    .ToArray();
                //for (int j = 0; j < xValuesBin.Length; j++)
                //{
                //    // check to see there is only one bin index in the x array
                    


                //    if (xValuesBin[binIndex] == null)
                //    {
                //        xValuesBin[binIndex] = new double[numSpectra];
                //        yValuesBin[binIndex] = new double[numSpectra];
                //    }

                //    xValuesBin[binIndex][i] = xArrays[i][j];
                //    yValuesBin[binIndex][i] = yArrays[i][j];
                //}
            }

            xValuesBin = xValuesBin.Where(p => p != null).ToArray();
            yValuesBin = yValuesBin.Where(p => p != null).ToArray();

            double[] xArray = new double[xValuesBin.Length];


            for (var i = 0; i < xArray.Length; i++)
            {
                xArray[i] = xValuesBin[i].Where(p => p != 0).Average();
                double x = xArray[i];
                var pixelStack = new PixelStack(x);
                pixelStack.AddIntensityVals(yValuesBin[i]);
                PixelStacks.Add(pixelStack);
            }
        }

        public void PerformNormalization(double[][] yArrays, double[] tics)
        {
            for (int i = 0; i < yArrays.Length; i++)
            {
                SpectrumNormalization.NormalizeSpectrumToTic(yArrays[i],
                    tics[i], tics.Average());
            }
        }

        public void CalculateNoiseEstimates(double[][] yArrays)
        {
            for (var i = 0; i < yArrays.Length; i++)
            {
                double[] yArray = yArrays[i];
                bool success = MRSNoiseEstimator.MRSNoiseEstimation(yArray, 0.01, out double noiseEstimate);
                // if the MRS noise estimate fails to converge, go by regular standard deviation
                if (!success)
                {
                    noiseEstimate = BasicStatistics.CalculateStandardDeviation(yArray);
                }

                NoiseEstimates.Add(i, noiseEstimate);
            }
        }

        public void CalculateScaleEstimates()
        {
            // insert avgdev method code here
        }

        public void CalculateWeights()
        {
            foreach (var entry in NoiseEstimates)
            {
                var successScale = ScaleEstimates.TryGetValue(entry.Key,
                    out double scale);
                if (!successScale) continue;

                var successNoise = NoiseEstimates.TryGetValue(entry.Key,
                    out double noise);
                if (!successNoise) continue;

                double weight = 1d / Math.Pow((scale * noise), 2);

                Weights.Add(entry.Key, weight);
            }
        }

        public void MergeSpectra()
        {
            double[] weights = Weights.OrderBy(i => i.Key)
                .Select(i => i.Value)
                .ToArray();
            foreach (var pixelStack in PixelStacks)
            {
                pixelStack.Average(weights);
            }
        }

        public double[][] GetMergedSpectrum()
        {
            double[] xArray = PixelStacks.Select(i => i.Mz).ToArray();
            double[] yArray = PixelStacks.Select(i => i.MergedValue).ToArray();
            return new[] { xArray, yArray };
        }
    }
}
