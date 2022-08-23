using MassSpectrometry;

namespace SpectralAveraging
{
    public class MultiScanDataObject 
    {
        public double MinX { get; set; }
        public double MaxX { get; set; }
        public int ScansToProcess { get; set; }
        public double[][] XArrays { get; set; }
        public double[][] YArrays { get; set; }
        public double[] TotalIonCurrent { get; set; }
        public MzSpectrum CompositeSpectrum { get; set; }
        public double? AverageIonCurrent
        {
            get { return TotalIonCurrent.Average(); }
        }

        public MultiScanDataObject(List<SingleScanDataObject> scanList)
        {
            GetMinX(scanList);
            GetMaxX(scanList);
            ProcessDataList(scanList); 
        }

        private void ProcessDataList(List<SingleScanDataObject> scanList)
        {
            ScansToProcess = scanList.Count;
            XArrays = new double[ScansToProcess][];
            YArrays = new double[ScansToProcess][];
            TotalIonCurrent = new double[ScansToProcess];

            for (int i = 0; i < ScansToProcess; i++)
            {
                XArrays[i] = scanList[i].XArray;
                YArrays[i] = scanList[i].YArray;
                TotalIonCurrent[i] = scanList[i].TotalIonCurrent;
            }
        }

        private void GetMinX(List<SingleScanDataObject> scanList)
        {
            MinX = scanList.Select(i => i.MinX).Aggregate((c,d) => c < d ? c : d); 
        }
        private void GetMaxX(List<SingleScanDataObject> scanList)
        {
            MaxX = scanList.Select(i => i.MaxX).Aggregate((c, d) => c > d ? c : d);
        }
    }
}
