using System.ComponentModel;
using System.Runtime.CompilerServices;
using System.Xml.Serialization;
using SpectralAveraging.NoiseEstimates;

namespace SpectralAveraging;

public class PixelStack
{
    // binned spectrum -> pixel stack -> average combination
    // TODO: Mz should not be a single value unless you're printing to console or getting the result. 
    // TODO: Store only values in _pixels and create GetIntensities method so you're only storing in 
    // one place. 
    public int Length => _pixels.Count;
    public int NonRejectedLength => _pixels.Count(i => i.Rejected == false); 
    public double MergedIntensityValue { get; private set; }
    public double MergedMzValue { get; private set; }
    public List<double> Intensity => _pixels.Select(i => i.Intensity).ToList();
    public double MzAverage => _pixels.Where(i => i.Rejected == false).Average(i => i.Mz);
    public IEnumerable<double> UnrejectedIntensities => _pixels.Where(i => i.Rejected == false)
        .Select(i => i.Intensity);
    public IEnumerable<double> UnrejectedMzs => _pixels.Where(i => i.Rejected == false)
        .Select(i => i.Mz); 
    private List<Pixel> _pixels { get; set; }
    
    // Implement INotifyPropertyChanged -> listen for Intensity changing -> Call method to 
    // update _pixels. 
    
    public PixelStack(IEnumerable<double> xArray, IEnumerable<double> yArray)
    {
        _pixels = xArray.Zip(yArray, (mz, its) => (mz,its))
            .Select((m,n) => new Pixel(n, m.mz, m.its, rejected:false))
            .ToList();
        _pixels.Sort(new Pixel.PixelComparer());
    }

    public IEnumerable<double> GetIntensities()
    {
        return _pixels.Select(i => i.Intensity);
    }

    public IEnumerable<double> GetMzValues()
    {
        return _pixels.Select(i => i.Mz); 
    }

    public IEnumerable<double> GetUnrejectedIntValues()
    {
        return _pixels.Where(i => i.Rejected == false)
            .Select(i => i.Intensity); 
    }
    public IEnumerable<double> GetUnrejectedMzValues()
    {
        return _pixels.Where(i => i.Rejected == false)
            .Select(i => i.Mz);
    }

    public void Reject(int index)
    {
        _pixels[index].Rejected = true; 
    }

    public bool IsIndexRejected(int index)
    {
        return _pixels[index].Rejected;
    }

    public void ModifyPixelIntensity(int index, double value)
    {
        _pixels[index].Intensity = value; 
    }

    public void ModifyPixelMz(int index, double value)
    {
        _pixels[index].Intensity = value; 
    }

    public double GetIntensityAtIndex(int index)
    {
        return _pixels[index].Intensity; 
    }

    public double GetMzAtIndex(int index)
    {
        return _pixels[index].Mz; 
    }

    public double CalculateMzAverage()
    {
        return _pixels.Where(j => j.Rejected == false)
            .Average(j => j.Mz); 
    }
    public void PerformRejection(SpectralAveragingOptions options)
    {
        OutlierRejection.RejectOutliers(this, options);
    }

    public void Average(IDictionary<int, double> weightsDictionary)
    {
        double numerator = 0;
        double denominator = 0;

        _pixels.Sort(new Pixel.PixelComparer());

        foreach (var weight in weightsDictionary)
        {
            int index = _pixels.IndexOf(_pixels.Where(i => i.SpectraId == weight.Key).First()); 
            if (_pixels[index].Rejected == true) continue;

            numerator += weight.Value * _pixels[index].Intensity;
            denominator += weight.Value; 
        }
        MergedIntensityValue = numerator / denominator;
        MergedMzValue = _pixels.Select(i => i.Mz).Average(); 
    }
    
    // People using the code should only ever see the list of doubles that represent 
    // the intensity values, and the merged spectra. 
    // But the actual backing probably needs to be more complicated than that. 
    // Use a sorted dictionary as the store

    internal class PixelStackComparer: IComparer<PixelStack>
    {
        public int Compare(PixelStack x, PixelStack y)
        {
            return x.CalculateMzAverage().CompareTo(y.CalculateMzAverage());
        }
    }
}

public class Pixel
{
    public int SpectraId;
    public double Intensity;
    public double Mz; 
    public bool Rejected; 
    public Pixel(int spectraId, double mz, double intensity, bool rejected)
    {
        SpectraId = spectraId; 
        Intensity = intensity;
        Rejected = rejected;
        Mz = mz; 
    }

    public Pixel()
    {

    }

    internal class PixelComparer : IComparer<Pixel>
    {
        public int Compare(Pixel? x, Pixel? y)
        {
            return x.SpectraId.CompareTo(y.SpectraId); 
        }
    }
}
