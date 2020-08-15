using System;
using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class ContinuousWavelet : WaveletBase
    {
        private readonly Func<double, Complex> _func;

        public ContinuousWavelet(Func<double, Complex> func,
            string shortName, string familyName,
            double freq = 0d,
            double min = double.NegativeInfinity, double max = double.PositiveInfinity)
        : base (shortName, familyName)
        {
            _func = func ?? throw new ArgumentNullException(nameof(func));
            base.CentralFrequency = freq;
            Minimum = min;
            Maximum = max;
        }

        public virtual double Minimum { get; }

        public virtual double Maximum { get; }

        public virtual Complex Evaluate(double time) => _func(time);

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
