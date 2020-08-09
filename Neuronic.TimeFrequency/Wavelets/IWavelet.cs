using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    public interface IWavelet
    {
        Complex[] Evaluate(double min, double max, int count);
        void Evaluate(double min, double max, Complex[] values, int start, int count);
        Complex Energy { get; }
        double CentralFrequency { get; }
    }
}
