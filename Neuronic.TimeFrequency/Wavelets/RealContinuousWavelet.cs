using System;

namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Represents a wavelet defined by a continuous function f: R -> R.
    /// </summary>
    public class RealContinuousWavelet : ContinuousWavelet, IWavelet<double>
    {
        private readonly Func<double, double> _func;

        /// <summary>
        /// Initializes a new instance of the <see cref="RealContinuousWavelet"/> class.
        /// </summary>
        /// <param name="func">The function.</param>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The family name.</param>
        /// <param name="freq">The central frequency.</param>
        /// <param name="min">The minimum time.</param>
        /// <param name="max">The maximum ime.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="max"/> is less than <paramref name="min"/>.</exception>
        /// <exception cref="ArgumentNullException">func</exception>
        public RealContinuousWavelet(Func<double, double> func, 
            string shortName, string familyName, 
            double freq = 0, double min = Double.NegativeInfinity, double max = Double.PositiveInfinity) 
            : base(x => func(x), shortName, familyName, freq, min, max)
        {
            _func = func ?? throw new ArgumentNullException(nameof(func));
        }

        /// <summary>
        /// Evaluates the function at the specified time instant.
        /// </summary>
        /// <param name="time">The time.</param>
        /// <returns>The evaluated value.</returns>
        public new virtual double Evaluate(double time) => _func(time);

        /// <summary>
        /// Evaluates the wavelet function in all it's domain.
        /// </summary>
        /// <returns>
        /// The evaluated values.
        /// </returns>
        public new Signal<double> EvaluateDomain()
        {
            return Evaluate(Minimum, Maximum, 1024);
        }

        /// <summary>
        /// Evaluates the wavelet function in the specified range.
        /// </summary>
        /// <param name="min">The minimum time.</param>
        /// <param name="max">The maximum time.</param>
        /// <param name="count">The amount of samples.</param>
        /// <returns>
        /// The evaluated values.
        /// </returns>
        public new Signal<double> Evaluate(double min, double max, int count)
        {
            if (count <= 0 || max < min)
                return new Signal<double>(new double[0]);

            var values = new double[count];
            var signal = new Signal<double>(values, min, (count - 1) / (max - min));
            Evaluate(signal);
            return signal;
        }

        /// <summary>
        /// Evaluates the wavelet function in the range defined by the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        public void Evaluate(Signal<double> signal)
        {
            var step = signal.SamplingPeriod;
            var min = signal.Start;
            for (int i = 0; i < signal.Count; i++)
                signal[i] = Evaluate(min + i * step);
        }
    }
}