using System.Numerics;

namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Abstraction of a wavelet function of the specified type.
    /// </summary>
    /// <typeparam name="T">The type of the wavelet's values (complex or real).</typeparam>
    public interface IWavelet<T>
    {
        /// <summary>
        /// Evaluates the wavelet function in all it's domain.
        /// </summary>
        /// <returns>The evaluated values.</returns>
        Signal<T> EvaluateDomain();
        /// <summary>
        /// Evaluates the wavelet function in the specified range.
        /// </summary>
        /// <param name="min">The minimum time.</param>
        /// <param name="max">The maximum time.</param>
        /// <param name="count">The amount of samples.</param>
        /// <returns>The evaluated values.</returns>
        Signal<T> Evaluate(double min, double max, int count);
        /// <summary>
        /// Evaluates the wavelet function in the range defined by the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        void Evaluate(Signal<T> signal);
        /// <summary>
        /// Gets the wavelet's estimated central frequency.
        /// </summary>
        /// <remarks>
        /// The central frequency if the wavelet's predominant frequency band according to it's frequency spectrum.
        /// </remarks>
        double CentralFrequency { get; }
    }
}
