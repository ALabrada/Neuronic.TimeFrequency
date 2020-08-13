using System;
using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class ContinuousWavelet : WaveletBase
    {
        private readonly Func<double, Complex> _func;
        private double _centralFrequency;

        public ContinuousWavelet(Func<double, Complex> func,
            double freq = 0d, Complex? energy = null,
            double min = double.NegativeInfinity, double max = double.PositiveInfinity)
        {
            _func = func ?? throw new ArgumentNullException(nameof(func));
            _centralFrequency = freq;
            Minimum = min;
            Maximum = max;
        }

        public ContinuousWavelet(Func<double, double> func,
            double freq = 0d, double energy = 0d,
            double min = double.NegativeInfinity, double max = double.PositiveInfinity)
            : this(new Func<double, Complex>(x => func(x)), freq, energy, min, max)
        {
        }

        public override double CentralFrequency => _centralFrequency > 0 ? _centralFrequency : (_centralFrequency = EstimateCentralFrequency(Evaluate()));
        public virtual double Minimum { get; }

        public virtual double Maximum { get; }

        public Complex Evaluate(double time) => _func(time);

        public override void Evaluate(Signal<Complex> signal)
        {
            var step = signal.SamplingPeriod;
            var min = signal.Delay;
            for (int i = 0; i < signal.Count; i++)
                signal[i] = Evaluate(min + i * step);
        }

        public override Signal<Complex> Evaluate()
        {
            return Evaluate(Minimum, Maximum, 1024);
        }
    }
}
