using System;
using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Represents a wavelet defined by a continuous function.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Wavelets.WaveletBase" />
    public class ContinuousWavelet : WaveletBase
    {
        private readonly Func<double, Complex> _func;

        /// <summary>
        /// Initializes a new instance of the <see cref="ContinuousWavelet"/> class.
        /// </summary>
        /// <param name="func">The wavelet function.</param>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The family name.</param>
        /// <param name="freq">The central frequency.</param>
        /// <param name="min">The minimum time.</param>
        /// <param name="max">The maximum time.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="max"/> is less than <paramref name="min"/>.</exception>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="func"/> is <c>null</c>.</exception>
        public ContinuousWavelet(Func<double, Complex> func,
            string shortName, string familyName,
            double freq = 0d,
            double min = double.NegativeInfinity, double max = double.PositiveInfinity)
        : base (shortName, familyName)
        {
            if (max <= min) throw new ArgumentOutOfRangeException(nameof(max));
            _func = func ?? throw new ArgumentNullException(nameof(func));
            base.CentralFrequency = freq;
            Minimum = min;
            Maximum = max;
        }

        /// <summary>
        /// Gets the minimum time value.
        /// </summary>
        /// <remarks>
        /// The function is zero for all values in the range (-Inf,<see cref="Minimum"/>).
        /// </remarks>
        public virtual double Minimum { get; }

        /// <summary>
        /// Gets the maximum time value.
        /// </summary>
        /// <remarks>
        /// The function is zero for all values in the range (<see cref="Maximum"/>, Inf).
        /// </remarks>
        public virtual double Maximum { get; }

        /// <summary>
        /// Evaluates the wavelet function at the specified time instant.
        /// </summary>
        /// <param name="time">The time.</param>
        /// <returns>The evaluated value.</returns>
        public virtual Complex Evaluate(double time) => _func(time);

        /// <summary>
        /// Evaluates the wavelet function in the range defined by the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        public override void Evaluate(Signal<Complex> signal)
        {
            var step = signal.SamplingPeriod;
            var min = signal.Start;
            for (int i = 0; i < signal.Count; i++)
                signal[i] = Evaluate(min + i * step);
        }

        /// <summary>
        /// Evaluates the wavelet function in all it's domain.
        /// </summary>
        /// <returns>
        /// The evaluated values.
        /// </returns>
        public override Signal<Complex> EvaluateDomain()
        {
            return Evaluate(Minimum, Maximum, 1024);
        }
    }
}
