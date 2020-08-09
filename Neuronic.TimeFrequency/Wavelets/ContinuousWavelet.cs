using System;
using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class ContinuousWavelet : WaveletBase
    {   
        private readonly Func<double, Complex> _func;

        public ContinuousWavelet(Func<double, Complex> func, 
            double freq = 0d, Complex? energy = null, 
            double min = double.NegativeInfinity, double max = double.PositiveInfinity)
        {
            _func = func ?? throw new ArgumentNullException(nameof(func));
            CentralFrequency = freq;
            Minimum = min;
            Maximum = max;

            var count = 1024;
            if (double.IsInfinity(min))
                min = 0;
            if (double.IsInfinity(max))
                max = 10;
            Complex[] values = null;
            if (CentralFrequency <= 0)
            {                
                values = Evaluate(min, max, count);
                CentralFrequency = EstimateCentralFrequency(values, (max - min) / (count - 1));
            }
            if (energy.HasValue)
                Energy = energy.Value;
            else
            {
                values = values ?? Evaluate(min, max, count);
                Energy = EstimateEnergy(values, min, max);
            }
        }

        public ContinuousWavelet(Func<double, double> func, 
            double freq = 0d, double energy = 0d,
            double min = double.NegativeInfinity, double max = double.PositiveInfinity)
            : this (new Func<double, Complex>(x => func(x)), freq, energy, min, max)
        {
        }

        public override Complex Energy { get; }

        public override double CentralFrequency { get; }

        public virtual double Minimum { get; }

        public virtual double Maximum { get; }

        public Complex Evaluate(double time) => _func(time);

        public override void Evaluate(double min, double max, Complex[] values, int start, int count)
        {
            var step = (max - min) / (count - 1);
            for (int i = 0; i < count; i++)
                values[i + start] = Evaluate(min + i * step);
        }
    }
}
