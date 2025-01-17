﻿using System;
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
            if (ticVal == 0)
                return;

            for (int i = 0; i < intensityArray.Length; i++)
            {
                intensityArray[i] /= ticVal;
            }
        }

        public static void NormalizeSpectrumToTic(double[] intensityArray, double ticVal)
        {
            if (ticVal == 0)
                return;

            for (int i = 0; i < intensityArray.Length; i++)
            {
                intensityArray[i] /= ticVal;
            }
        }

        public static void NormalizeSpectrumToTic(double[] intensityArray, double ticVal, double avgTicVal)
        {
            if (ticVal == 0)
                return;

            for (int i = 0; i < intensityArray.Length; i++)
            {
                intensityArray[i] = intensityArray[i] / ticVal * avgTicVal;
            }
        }
    }
}
