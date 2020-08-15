namespace Neuronic.TimeFrequency.Wavelets
{
    public static class WaveletExtensions
    {
        public static double GetFrequencyOf(this IWavelet wavelet, double scale, double samplingPeriod = 1d)
        {
            return wavelet.CentralFrequency / (scale * samplingPeriod);
        }

        public static double GetScaleFor(this IWavelet wavelet, double frequency, double samplingPeriod = 1d)
        {
            return wavelet.CentralFrequency / (frequency * samplingPeriod);
        }
    }
}
