using System;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class RealContinuousWavelet : ContinuousWavelet, IWavelet<double>
    {
        private readonly Func<double, double> _func;

        public RealContinuousWavelet(Func<double, double> func, 
            string shortName, string familyName, 
            double freq = 0, double min = Double.NegativeInfinity, double max = Double.PositiveInfinity) 
            : base(x => func(x), shortName, familyName, freq, min, max)
        {
            _func = func;
        }

        public new virtual double Evaluate(double time) => _func(time);

        public new Signal<double> Evaluate()
        {
            return Evaluate(Minimum, Maximum, 1024);
        }

        public new Signal<double> Evaluate(double min, double max, int count)
        {
            var values = new double[count];
            var signal = new Signal<double>(values, min, (count - 1) / (max - min));
            Evaluate(signal);
            return signal;
        }

        public void Evaluate(Signal<double> signal)
        {
            var step = signal.SamplingPeriod;
            var min = signal.Delay;
            for (int i = 0; i < signal.Count; i++)
                signal[i] = Evaluate(min + i * step);
        }
    }
}