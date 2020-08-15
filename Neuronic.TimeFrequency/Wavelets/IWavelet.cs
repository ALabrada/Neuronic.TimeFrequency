using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    public interface IWavelet
    {
        Signal<Complex> Evaluate();
        Signal<Complex> Evaluate(double min, double max, int count);
        void Evaluate(Signal<Complex> signal);
        double CentralFrequency { get; }
    }
}
