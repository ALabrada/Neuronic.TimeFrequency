using System;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;

namespace Neuronic.TimeFrequency.Wavelets
{
    public abstract class WaveletBase : IWavelet
    {
        private Complex _energy;
        private double _centralFrequency;

        public virtual Complex Energy
        {
            get => _energy != Complex.Zero ? _energy : (_energy = EstimateEnergy(Evaluate()));
            protected set => _energy = value;
        }

        protected Complex EstimateEnergy(Complex[] values, double min, double step)
        {
            var sum = Complex.Zero;
            for (int i = 0; i < values.Length; i++)
            {
                var x = min + step * i;
                var y = values[i];
                sum += step * y * y / x;
            }
            return 2 * Math.PI * sum;
        }

        protected Complex EstimateEnergy(Signal<Complex> signal)
        {
            return EstimateEnergy(signal.Samples, signal.Delay, signal.SamplingPeriod);
        }

        protected double EstimateEnergy(double[] values, double min, double step)
        {
            var sum = 0d;
            for (int i = 0; i < values.Length; i++)
            {
                var x = min + step * i;
                var y = values[i];
                sum += step * y * y / x;
            }
            return 2 * Math.PI * sum;
        }

        protected double EstimateEnergy(Signal<double> signal)
        {
            return EstimateEnergy(signal.Samples, signal.Delay, signal.SamplingPeriod);
        }

        public virtual double CentralFrequency
        {
            get => _centralFrequency >= 0 ? _centralFrequency : (_centralFrequency = EstimateCentralFrequency(Evaluate()));
            protected set => _centralFrequency = value;
        }

        protected double EstimateCentralFrequency(Complex[] values, double samplingPeriod)
        {
            FourierTransform2.FFT(values, Accord.Math.FourierTransform.Direction.Forward);
            var maxIndex = 0;
            for (int i = 1; i < values.Length / 2; i++)
                if (values[i].Magnitude > values[maxIndex].Magnitude)
                    maxIndex = i;
            var domain = values.Length * samplingPeriod;
            return maxIndex / domain;
        }

        protected double EstimateCentralFrequency(Signal<Complex> signal) => EstimateCentralFrequency(signal.Samples, signal.SamplingPeriod);

        protected double EstimateCentralFrequency(double[] values, double samplingPeriod)
        {
            return EstimateCentralFrequency(values.ToComplex(), samplingPeriod);
        }

        protected double EstimateCentralFrequency(Signal<double> signal) => EstimateCentralFrequency(signal.Samples, signal.SamplingPeriod);

        public virtual Signal<Complex> Evaluate(double min, double max, int count)
        {
            var values = new Complex[count];
            var signal = new Signal<Complex>(values, min, (count - 1) / (max - min));
            Evaluate(signal);
            return signal;
        }

        public abstract void Evaluate(Signal<Complex> signal);
        public abstract Signal<Complex> Evaluate();
    }
}
