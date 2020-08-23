using System;
using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Base class for wavelet functions.  
    /// </summary>
    public abstract class WaveletBase : IWavelet<Complex>
    {
        private double _centralFrequency;

        /// <summary>
        /// Initializes a new instance of the <see cref="WaveletBase"/> class.
        /// </summary>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The wavelet family name.</param>
        protected WaveletBase(string shortName, string familyName)
        {
            ShortName = shortName;
            FamilyName = familyName;
        }

        /// <summary>
        /// Gets the short name.
        /// </summary>
        public string ShortName { get; }

        /// <summary>
        /// Gets the name of the wavelet family.
        /// </summary>
        public string FamilyName { get; }

        /// <summary>
        /// Gets the wavelet's estimated central frequency.
        /// </summary>
        /// <remarks>
        /// The central frequency if the wavelet's predominant frequency band according to it's frequency spectrum.
        /// </remarks>
        public virtual double CentralFrequency
        {
            get => _centralFrequency > 0 ? _centralFrequency : (_centralFrequency = EstimateCentralFrequency(EvaluateDomain()));
            protected set => _centralFrequency = value;
        }

        /// <summary>
        /// Estimates the central frequency from the specified values.
        /// </summary>
        /// <param name="values">The values.</param>
        /// <param name="samplingPeriod">The sampling period.</param>
        /// <returns>The estimated central frequency.</returns>
        protected double EstimateCentralFrequency(Complex[] values, double samplingPeriod)
        {
            values.FFT();
            var maxIndex = 0;
            for (int i = 1; i < values.Length / 2; i++)
                if (values[i].Magnitude > values[maxIndex].Magnitude)
                    maxIndex = i;
            var domain = values.Length * samplingPeriod;
            return maxIndex / domain;
        }

        /// <summary>
        /// Estimates the central frequency from the specified values.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <returns>The estimated central frequency.</returns>
        protected double EstimateCentralFrequency(Signal<Complex> signal) => EstimateCentralFrequency(signal.Samples, signal.SamplingPeriod);

        /// <summary>
        /// Estimates the central frequency from the specified values.
        /// </summary>
        /// <param name="values">The values.</param>
        /// <param name="samplingPeriod">The sampling period.</param>
        /// <returns>The estimated central frequency.</returns>
        protected double EstimateCentralFrequency(double[] values, double samplingPeriod)
        {
            return EstimateCentralFrequency(values.ToComplex(), samplingPeriod);
        }

        /// <summary>
        /// Estimates the central frequency from the specified values.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <returns>The estimated central frequency.</returns>
        protected double EstimateCentralFrequency(Signal<double> signal) => EstimateCentralFrequency(signal.Samples, signal.SamplingPeriod);

        /// <summary>
        /// Evaluates the wavelet function in the specified range.
        /// </summary>
        /// <param name="min">The minimum time.</param>
        /// <param name="max">The maximum time.</param>
        /// <param name="count">The amount of samples.</param>
        /// <returns>
        /// The evaluated values.
        /// </returns>
        public virtual Signal<Complex> Evaluate(double min, double max, int count)
        {
            if (count <= 0 || max < min)
                return new Signal<Complex>(new Complex[0]);

            var values = new Complex[count];
            var signal = new Signal<Complex>(values, min, (count - 1) / (max - min));
            Evaluate(signal);
            return signal;
        }

        /// <summary>
        /// Evaluates the wavelet function in the range defined by the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        public abstract void Evaluate(Signal<Complex> signal);

        /// <summary>
        /// Evaluates the wavelet function in all it's domain.
        /// </summary>
        /// <returns>
        /// The evaluated values.
        /// </returns>
        public abstract Signal<Complex> EvaluateDomain();
    }
}
