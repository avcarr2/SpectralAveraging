using MassSpectrometry;

namespace SpectralAveraging
{
    public class SingleScanDataObject 
    {
        public double[] XArray { get; private set; }
        public double[] YArray { get; private set; }
        public double TotalIonCurrent { get; private set; }
        public double MinX { get; private set; }
        public double MaxX { get; private set; }

        public double Resolution { get; private set; }


        public SingleScanDataObject(MsDataScan scan)
        {
            XArray = scan.MassSpectrum.XArray;
            YArray = scan.MassSpectrum.YArray;
            TotalIonCurrent = scan.TotalIonCurrent;
            MinX = scan.MassSpectrum.XArray.Min();
            MaxX = scan.MassSpectrum.XArray.Max();
        }
        public void UpdateYarray(double[] newYarray)
        {
            YArray = newYarray; 
        }

        // TODO: Target for optimization 
        public static List<SingleScanDataObject> ConvertMSDataScansInBulk(List<MsDataScan> scans)
        {
            List<SingleScanDataObject> singleScanDataObjects = new List<SingleScanDataObject>();
            foreach (var scan in scans)
            {
                singleScanDataObjects.Add(new SingleScanDataObject(scan));
            }
            return singleScanDataObjects;
        }

        // temp
        public List<double>[] GetChargeStateEnvelopeMz()
        {
            List<double>[] result = new List<double>[1] { new List<double>()};
            return result;
        }
    }
}
