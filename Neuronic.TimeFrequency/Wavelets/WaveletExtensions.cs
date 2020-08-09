namespace Neuronic.TimeFrequency.Wavelets
{
    public static class WaveletExtensions
    {
        public static double GetFrequencyOf(this IWavelet wavelet, double scale)
        {
            return wavelet.CentralFrequency / scale;
        }

        public static double GetScaleFor(this IWavelet wavelet, double frequency)
        {
            return wavelet.CentralFrequency / frequency;
        }
    }
}
