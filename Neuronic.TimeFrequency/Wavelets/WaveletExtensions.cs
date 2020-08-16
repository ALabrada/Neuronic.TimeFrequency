namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Contains extension methods for the wavelets.
    /// </summary>
    public static class WaveletExtensions
    {
        /// <summary>
        /// Gets the predominant frequency of the wavelet at the specified scale.
        /// </summary>
        /// <typeparam name="T">The type of the wavelet values.</typeparam>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="scale">The scale.</param>
        /// <param name="samplingPeriod">The sampling period.</param>
        /// <returns>The predominant frequency.</returns>
        public static double GetFrequencyOf<T>(this IWavelet<T> wavelet, double scale, double samplingPeriod = 1d)
        {
            return wavelet.CentralFrequency / (scale * samplingPeriod);
        }

        /// <summary>
        /// Gets the best scale for scaling the wavelet to the specified frequency..
        /// </summary>
        /// <typeparam name="T">The type of the wavelet values.</typeparam>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="frequency">The frequency.</param>
        /// <param name="samplingPeriod">The sampling period.</param>
        /// <returns>The predominant frequency.</returns>
        public static double GetScaleFor<T>(this IWavelet<T> wavelet, double frequency, double samplingPeriod = 1d)
        {
            return wavelet.CentralFrequency / (frequency * samplingPeriod);
        }
    }
}
