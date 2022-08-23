using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectralAveraging
{
    public static class SpectrumNormalization
    {
        public static void NormalizeSpectrumToTic(ref double[] intensityArray, double ticVal)
        {
            for (int i = 0; i < intensityArray.Length; i++)
            {
                intensityArray[i] /= ticVal;
            }
        }

        public static void NormalizeSpectrumToTic(double[] intensityArray, double ticVal)
        {
            for (int i = 0; i < intensityArray.Length; i++)
            {
                intensityArray[i] /= ticVal;
            }
        }

        public static void NormalizeSpectrumToTic(double[] intensityArray, double ticVal, double avgTicVal)
        {
            for (int i = 0; i < intensityArray.Length; i++)
            {
                intensityArray[i] = intensityArray[i] / ticVal * avgTicVal;
            }
        }

        /// <summary>
        /// Normalization method overload for MultiScanDataObjects
        /// </summary>
        /// <param name="scans"></param>
        public static void NormalizeSpectrumToTic(MultiScanDataObject scans, bool multiplyByAverageTic)
        {
            for (int i = 0; i < scans.YArrays.GetLength(0); i++)
            {
                if (multiplyByAverageTic && scans.AverageIonCurrent != null)
                {
                    NormalizeSpectrumToTic(scans.YArrays[i], scans.TotalIonCurrent[i], (double)scans.AverageIonCurrent);
                }
                else
                {
                    NormalizeSpectrumToTic(scans.YArrays[i], scans.TotalIonCurrent[i]);
                }
            }
        }

        /// <summary>
        /// Normalization method overload for SingleScanDataObjects
        /// </summary>
        /// <param name="scan"></param>
        public static void NormalizeSpectrumToTic(SingleScanDataObject scan)
        {
            NormalizeSpectrumToTic(scan.YArray, scan.TotalIonCurrent);
        }
    }
}
