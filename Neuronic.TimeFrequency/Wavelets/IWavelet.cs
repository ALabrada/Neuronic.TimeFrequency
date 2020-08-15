using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    public interface IWavelet<T>
    {
        Signal<T> Evaluate();
        Signal<T> Evaluate(double min, double max, int count);
        void Evaluate(Signal<T> signal);
        double CentralFrequency { get; }
    }
}
