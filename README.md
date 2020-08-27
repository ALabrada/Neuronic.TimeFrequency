# Neuronic.TimeFrequency
Feature extraction library specialized in Time-Frequency domain features. Includes Wavelet transforms, spectrograms, time-frequency distributions and other algorithms.

## Code Examples

### Spectrogram using Short-Time Fourier Transform
```
var window = MathNet.Numerics.Window.Gauss(51, 0.5);
var spectrogram = Spectrogram.Estimate(new Signal<double>(samples, fs: 1000d), window);
Console.WriteLine("Amplitude: {0}", spectrogram.EnumerateValueOfFrequency(60d).Max());
```

### Continuous Wavelet Transform
```
var wavelet = Wavelets.Wavelets.Morlet;
var scales = Enumerable.Range(0, 10);
var cwt = ContinuousWaveletTransform.EstimateUsingConvolutions(new Signal<double>(samples, fs: 1000d), wavelet, scales);
Console.WriteLine("Amplitude: {0}", cwt.EnumerateValuesOfFrequencyAt(4).Max());
```

### Choi-Williams Distribution
```
var kernel = new ChoiWilliamsDistribution {Sigma = 30};
var cwd = TimeFrequencyDistribution.Estimate(new Signal<double>(samples, fs: 1000d), kernel);
Console.WriteLine("Amplitude: {0}", cwd.EnumerateValueOfFrequency(60d).Max());
```

### Hilbert-Huang Transform
```
var emd = EmpiricalModeDecomposition.Estimate(new Signal<double>(samples, fs: 1000d));
var hht = emd.HilbertSpectralAnalysis();
Console.WriteLine("Frequency: {0}", hht[0].Frequency.Min());
```

## Instalation
The binaries are available at [Nuget](https://www.nuget.org/packages/Neuronic.TimeFrequency/). To install it run the command `Install-Package Neuronic.TimeFrequency`.

## API Reference
The previously mentioned algorithms are located in the `Neuronic.TimeFrequency.Transforms` namespace and
implement the `ITimeFrequencyRepresentation` interface. It's main feature is and indexer that obtains the
amplitude for a given coordenate in the time and frequency axis. Most of them also implement the more 
specialized IBilinearTimeFrequencyRepresentation, that exposes the inner sampling of the frequency domain.

The `ContinuousWaveletTransform` can employ two different algorithms: the first based on the Fast Fourier Transform
algorithm [1] and the second on the convolution-based implementation from MATLAB. Both can use continuous
or orthogonal wavelets, available in the `Neuronic.TimeFrequency.Wavelets` namespace. The `Wavelets` static class
exposes some of the most used ones.

The `TimeFrequencyDistribution` algorithm is based on [2] and can use differnte times of kernels, 
including the `ChoiWilliamsDistribution`. The predefined kernels are available in the 
`Neuronic.TimeFrequency.Kernels` namespace.

The `EmpiricalModeDecomposition` implementation is based on [3]. It also includes an implementation of 
the Hilbert Spectral Analysis and the spectral analysis algorithm proposed in [3].

## References
1. Stark, Hans-Georg: "Wavelets and Signal Processing: An Application-Based Introduction", Springer, 2005.
2. Toole, J.M., Boashash, B.: "Fast and memory-efficient algorithms for computing quadratic time–frequency distributions", Applied and Computational Harmonic Analysis, vol. 35, no. 2, pp. 350–358, 2013.
3. Rato, R. T., Ortigueira, M. D. and Batista, A. G.: "On the HHT, its problems, and some solutions", Mechanical Systems and Signal Processing , vol. 22, no. 6, pp. 1374-1394, August 2008.