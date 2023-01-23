using Plotly.NET;

namespace Tests;


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
    public static (double[][] xArrays, double[][] yArrays) GetJaggedArrays(this List<SimulatedProtein> peakList)
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
    public static (double[][] xArrays, double[][] yArrays) GetJaggedArrays(this List<TestData> peakList)
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
