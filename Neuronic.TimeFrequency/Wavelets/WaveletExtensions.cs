namespace Neuronic.TimeFrequency.Wavelets
{
    public static class WaveletExtensions
    {
        public static double GetFrequencyOf<T>(this IWavelet<T> wavelet, double scale, double samplingPeriod = 1d)
        {
            return wavelet.CentralFrequency / (scale * samplingPeriod);
        }

        public static double GetScaleFor<T>(this IWavelet<T> wavelet, double frequency, double samplingPeriod = 1d)
        {
            return wavelet.CentralFrequency / (frequency * samplingPeriod);
        }
    }
}
