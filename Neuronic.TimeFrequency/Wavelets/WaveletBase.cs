using System;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;

namespace Neuronic.TimeFrequency.Wavelets
{
    public abstract class WaveletBase : IWavelet
    {
        public abstract Complex Energy { get; }

        protected Complex EstimateEnergy(Complex[] values, double min, double max)
        {
            var step = (max - min) / (values.Length - 1);
            var sum = Complex.Zero;
            for (int i = 0; i < values.Length; i++)
            {
                var x = min + step * i;
                var y = values[i];
                sum += step * y * y / x;
            }
            return 2 * Math.PI * sum;
        }

        protected Complex EstimateEnergy(double[] values, double min, double max)
        {
            var step = (max - min) / (values.Length - 1);
            var sum = 0d;
            for (int i = 0; i < values.Length; i++)
            {
                var x = min + step * i;
                var y = values[i];
                sum += step * y * y / x;
            }
            return 2 * Math.PI * sum;
        }

        public abstract double CentralFrequency { get; }

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

        protected double EstimateCentralFrequency(double[] values, double samplingPeriod)
        {
            return EstimateCentralFrequency(values.ToComplex(), samplingPeriod);
        }

        public Complex[] Evaluate(double min, double max, int count)
        {
            var result = new Complex[count];
            Evaluate(min, max, result, 0, count);
            return result;
        }

        public abstract void Evaluate(double min, double max, Complex[] values, int start, int count);
    }
}
