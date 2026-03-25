const palette = {
  sc: "#0f766e",
  mrc: "#b45309",
  sskmSc: "#7c3aed",
  sskmMrc: "#b42318",
  siso: "#0f766e",
  ref: "#0f766e",
  rx: "#b45309",
};

const presets = {
  assignment2: {
    fftSize: 128,
    giRatio: 0.25,
    modulationOrder: 4,
    multiPath: 7,
    snrStart: 0,
    snrEnd: 30,
    snrStep: 3,
  },
  assignment3: {
    fftSize: 128,
    giRatio: 0.25,
    modulationOrder: 4,
    multiPath: 7,
    snrStart: 0,
    snrEnd: 30,
    snrStep: 3,
  },
};

const dom = {
  scenario: document.getElementById("scenario"),
  berTargetIterations: document.getElementById("berTargetIterations"),
  berTargetIterationsValue: document.getElementById("berTargetIterationsValue"),
  seed: document.getElementById("seed"),
  runButton: document.getElementById("runButton"),
  presetButton: document.getElementById("presetButton"),
  berTopControls: document.getElementById("berTopControls"),
  chartCanvas: document.getElementById("chartCanvas"),
  constellationCanvas: document.getElementById("constellationCanvas"),
  legend: document.getElementById("legend"),
  liveCards: document.getElementById("liveCards"),
  fixedConfig: document.getElementById("fixedConfig"),
  resultsTable: document.getElementById("resultsTable"),
  resultsHeader: document.getElementById("resultsHeader"),
  statusText: document.getElementById("statusText"),
  assignment1Status: document.getElementById("assignment1Status"),
  assignment1Metrics: document.getElementById("assignment1Metrics"),
  berSection: document.getElementById("berSection"),
  assignment1Section: document.getElementById("assignment1Section"),
  a1SampleCount: document.getElementById("a1SampleCount"),
  a1ChannelGain: document.getElementById("a1ChannelGain"),
  a1ChannelPhase: document.getElementById("a1ChannelPhase"),
  a1NoiseSigma: document.getElementById("a1NoiseSigma"),
  a1SampleCountValue: document.getElementById("a1SampleCountValue"),
  a1ChannelGainValue: document.getElementById("a1ChannelGainValue"),
  a1ChannelPhaseValue: document.getElementById("a1ChannelPhaseValue"),
  a1NoiseSigmaValue: document.getElementById("a1NoiseSigmaValue"),
};

let runToken = 0;

class RNG {
  constructor(seed) {
    this.state = (seed >>> 0) || 1;
  }

  next() {
    this.state = (1664525 * this.state + 1013904223) >>> 0;
    return this.state / 4294967296;
  }

  normalPair() {
    const u1 = Math.max(this.next(), 1e-12);
    const u2 = this.next();
    const mag = Math.sqrt(-2 * Math.log(u1));
    const angle = 2 * Math.PI * u2;
    return [mag * Math.cos(angle), mag * Math.sin(angle)];
  }
}

function complex(re = 0, im = 0) {
  return { re, im };
}

function add(a, b) {
  return { re: a.re + b.re, im: a.im + b.im };
}

function sub(a, b) {
  return { re: a.re - b.re, im: a.im - b.im };
}

function mul(a, b) {
  return {
    re: a.re * b.re - a.im * b.im,
    im: a.re * b.im + a.im * b.re,
  };
}

function scale(a, value) {
  return { re: a.re * value, im: a.im * value };
}

function conj(a) {
  return { re: a.re, im: -a.im };
}

function abs2(a) {
  return a.re * a.re + a.im * a.im;
}

function div(a, b) {
  const denom = abs2(b) + 1e-12;
  return {
    re: (a.re * b.re + a.im * b.im) / denom,
    im: (a.im * b.re - a.re * b.im) / denom,
  };
}

function fft(input) {
  const n = input.length;
  if (n <= 1) return input.slice();

  const even = [];
  const odd = [];
  for (let i = 0; i < n; i += 1) {
    (i % 2 === 0 ? even : odd).push(input[i]);
  }

  const fftEven = fft(even);
  const fftOdd = fft(odd);
  const result = new Array(n);

  for (let k = 0; k < n / 2; k += 1) {
    const angle = (-2 * Math.PI * k) / n;
    const twiddle = complex(Math.cos(angle), Math.sin(angle));
    const t = mul(twiddle, fftOdd[k]);
    result[k] = add(fftEven[k], t);
    result[k + n / 2] = sub(fftEven[k], t);
  }

  return result;
}

function ifft(input) {
  return fft(input.map(conj)).map(conj).map((value) => scale(value, 1 / input.length));
}

function convolve(signal, channel) {
  const out = Array.from({ length: signal.length + channel.length - 1 }, () => complex(0, 0));
  for (let i = 0; i < signal.length; i += 1) {
    for (let j = 0; j < channel.length; j += 1) {
      out[i + j] = add(out[i + j], mul(signal[i], channel[j]));
    }
  }
  return out;
}

function createSnrRange(start, end, step) {
  const values = [];
  for (let snr = start; snr <= end + 1e-9; snr += step) values.push(Number(snr.toFixed(6)));
  return values;
}

function randomBits(length, rng) {
  return Array.from({ length }, () => (rng.next() > 0.5 ? 1 : 0));
}

function rayleighChannel(multiPath, rng) {
  const profile = [];
  let total = 0;
  for (let i = 0; i < multiPath; i += 1) {
    const power = Math.exp(-(i + 1) / 5);
    profile.push(power);
    total += power;
  }

  return profile.map((power) => {
    const [n1, n2] = rng.normalPair();
    const sigma = Math.sqrt(0.5 * (power / total));
    return complex(n1 * sigma, n2 * sigma);
  });
}

function addAwgn(signal, snrDb, rng) {
  const noisePower = 10 ** (-snrDb / 10);
  const sigma = Math.sqrt(noisePower / 2);
  return signal.map((sample) => {
    const [n1, n2] = rng.normalPair();
    return add(sample, complex(n1 * sigma, n2 * sigma));
  });
}

function qpskMod(bits) {
  const symbols = [];
  for (let i = 0; i < bits.length; i += 2) {
    const imag = bits[i] ? 1 : -1;
    const real = bits[i + 1] ? 1 : -1;
    symbols.push(complex(real * 0.7071, imag * 0.7071));
  }
  return symbols;
}

function qpskDemod(symbols) {
  const bits = [];
  symbols.forEach((symbol) => {
    bits.push(symbol.im > 0 ? 1 : 0);
    bits.push(symbol.re > 0 ? 1 : 0);
  });
  return bits;
}

function qam16Mod(bits) {
  const symbols = [];
  for (let i = 0; i < bits.length; i += 4) {
    const b1 = bits[i];
    const b2 = bits[i + 1];
    const b3 = bits[i + 2];
    const b4 = bits[i + 3];
    const first = b1 * 4 - 2;
    const second = (b2 !== b1 ? 1 : 0) * 2 - 1;
    const third = b3 * 4 - 2;
    const fourth = (b4 !== b3 ? 1 : 0) * 2 - 1;
    symbols.push(complex((first + second) * 0.3162, (third + fourth) * 0.3162));
  }
  return symbols;
}

function qam16Demod(symbols) {
  const bits = [];
  symbols.forEach((symbol) => {
    bits.push(symbol.re > 0 ? 1 : 0);
    bits.push(Math.abs(symbol.re) < 0.6325 ? 1 : 0);
    bits.push(symbol.im > 0 ? 1 : 0);
    bits.push(Math.abs(symbol.im) < 0.6325 ? 1 : 0);
  });
  return bits;
}

function baseMod(bits, modulationOrder) {
  return modulationOrder === 2 ? qpskMod(bits) : qam16Mod(bits);
}

function baseDemod(symbols, modulationOrder) {
  return modulationOrder === 2 ? qpskDemod(symbols) : qam16Demod(symbols);
}

function nearestConstellationDemod(symbols, modulationOrder) {
  const bitsPerSymbol = modulationOrder;
  const referenceBits = [];
  const totalSymbols = 2 ** bitsPerSymbol;

  for (let value = 0; value < totalSymbols; value += 1) {
    const symbolBits = [];
    for (let bitIndex = bitsPerSymbol - 1; bitIndex >= 0; bitIndex -= 1) {
      symbolBits.push((value >> bitIndex) & 1);
    }
    referenceBits.push(symbolBits);
  }

  const flattened = referenceBits.flat();
  const references = baseMod(flattened, modulationOrder);
  const demodBits = [];

  symbols.forEach((symbol) => {
    let bestIndex = 0;
    let bestDistance = Number.POSITIVE_INFINITY;
    for (let index = 0; index < references.length; index += 1) {
      const diffRe = symbol.re - references[index].re;
      const diffIm = symbol.im - references[index].im;
      const distance = diffRe * diffRe + diffIm * diffIm;
      if (distance < bestDistance) {
        bestDistance = distance;
        bestIndex = index;
      }
    }
    demodBits.push(...referenceBits[bestIndex]);
  });

  return demodBits;
}

function addCyclicPrefix(ofdm, giSize) {
  return ofdm.slice(ofdm.length - giSize).concat(ofdm);
}

function removeCyclicPrefix(rx, giSize, fftSize) {
  return rx.slice(giSize, giSize + fftSize);
}

function normalizeIfft(symbols) {
  return ifft(symbols).map((value) => scale(value, Math.sqrt(symbols.length)));
}

function normalizeFft(samples) {
  return fft(samples).map((value) => scale(value, 1 / Math.sqrt(samples.length)));
}

function channelResponse(channel, fftSize) {
  const padded = Array.from({ length: fftSize }, (_, index) => channel[index] || complex(0, 0));
  return fft(padded);
}

function countBitErrors(a, b) {
  let count = 0;
  for (let i = 0; i < a.length; i += 1) {
    if (a[i] !== b[i]) count += 1;
  }
  return count;
}

function equalize(received, channel) {
  return received.map((value, index) => div(value, channel[index]));
}

function selectCombine(receivedA, receivedB, channelA, channelB) {
  const y = [];
  const h = [];
  for (let i = 0; i < receivedA.length; i += 1) {
    if (abs2(channelA[i]) >= abs2(channelB[i])) {
      y.push(receivedA[i]);
      h.push(channelA[i]);
    } else {
      y.push(receivedB[i]);
      h.push(channelB[i]);
    }
  }
  return { y, h };
}

function mrcCombine(receivedA, receivedB, channelA, channelB) {
  const combined = [];
  const denom = [];
  for (let i = 0; i < receivedA.length; i += 1) {
    combined.push(add(mul(receivedA[i], conj(channelA[i])), mul(receivedB[i], conj(channelB[i]))));
    denom.push(complex(abs2(channelA[i]) + abs2(channelB[i]), 0));
  }
  return { y: combined, h: denom };
}

function encodeSskm(bits) {
  const encoded = new Array(bits.length * 2);
  for (let i = 0; i < bits.length; i += 1) {
    encoded[2 * i] = bits[i];
    encoded[2 * i + 1] = bits[i] ? 0 : 1;
  }
  return encoded;
}

function decodeSskm(bits) {
  const decoded = new Array(bits.length / 2);
  for (let i = 0; i < decoded.length; i += 1) {
    decoded[i] = bits[2 * i] > bits[2 * i + 1] ? 1 : 0;
  }
  return decoded;
}

function getAssignment1Config() {
  return {
    modulationOrder: 4,
    sampleCount: Number(dom.a1SampleCount.value),
    channelGain: Number(dom.a1ChannelGain.value),
    channelPhaseDeg: Number(dom.a1ChannelPhase.value),
    noiseSigma: Number(dom.a1NoiseSigma.value),
    seed: Number(dom.seed.value),
  };
}

function getBerConfig() {
  const scenario = dom.scenario.value;
  return {
    scenario,
    ...presets[scenario],
    targetIterations: Number(dom.berTargetIterations.value),
    seed: Number(dom.seed.value),
  };
}

function updateAssignment1Labels() {
  dom.a1SampleCountValue.textContent = dom.a1SampleCount.value;
  dom.a1ChannelGainValue.textContent = Number(dom.a1ChannelGain.value).toFixed(2);
  dom.a1ChannelPhaseValue.textContent = dom.a1ChannelPhase.value;
  dom.a1NoiseSigmaValue.textContent = Number(dom.a1NoiseSigma.value).toFixed(2);
}

function updateBerTargetLabel() {
  dom.berTargetIterationsValue.textContent = dom.berTargetIterations.value;
}

function simulateAssignment1() {
  const config = getAssignment1Config();
  const rng = new RNG(config.seed);
  const idealBits = [];
  const symbolTemplates = config.modulationOrder === 2
    ? [
        [0, 0],
        [0, 1],
        [1, 0],
        [1, 1],
      ]
    : [
        [0, 0, 0, 0],
        [0, 1, 0, 1],
        [1, 1, 1, 1],
        [1, 0, 1, 0],
        [0, 0, 1, 1],
        [0, 1, 1, 0],
        [1, 1, 0, 0],
        [1, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 1, 0, 0],
        [1, 1, 1, 0],
        [1, 0, 1, 1],
        [0, 0, 1, 0],
        [0, 1, 1, 1],
        [1, 1, 0, 1],
        [1, 0, 0, 0],
      ];

  symbolTemplates.forEach((entry) => idealBits.push(...entry));
  const references = baseMod(idealBits, config.modulationOrder);

  const randomBitsForRx = randomBits(Math.max(config.modulationOrder, config.sampleCount), rng);
  const transmitted = baseMod(randomBitsForRx, config.modulationOrder);
  const phaseRad = (config.channelPhaseDeg * Math.PI) / 180;
  const channel = complex(config.channelGain * Math.cos(phaseRad), config.channelGain * Math.sin(phaseRad));
  const received = transmitted.map((symbol) => {
    const faded = mul(channel, symbol);
    const [n1, n2] = rng.normalPair();
    return add(faded, complex(n1 * config.noiseSigma, n2 * config.noiseSigma));
  });

  const mappedBits = nearestConstellationDemod(received, config.modulationOrder).slice(0, randomBitsForRx.length);
  const bitErrors = countBitErrors(randomBitsForRx, mappedBits);
  const ber = bitErrors / randomBitsForRx.length;

  drawConstellation(references, received);
  dom.assignment1Status.textContent = `H = ${channel.re.toFixed(2)} ${channel.im >= 0 ? "+" : "-"} ${Math.abs(channel.im).toFixed(2)}j`;
  dom.assignment1Metrics.innerHTML = `
    <div class="summary-card">
      <div>Mapped BER</div>
      <strong>${ber < 0.001 ? ber.toExponential(2) : ber.toFixed(4)}</strong>
      <div>Nearest reference point demodulation</div>
    </div>
    <div class="summary-card">
      <div>Bit Errors</div>
      <strong>${bitErrors}</strong>
      <div>Out of ${randomBitsForRx.length} bits</div>
    </div>
    <div class="summary-card">
      <div>Reference Size</div>
      <strong>${references.length}</strong>
      <div>Fixed ideal constellation points</div>
    </div>
  `;
}

function drawConstellation(references, received) {
  const canvas = dom.constellationCanvas;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  const pad = 48;
  const points = references.concat(received);
  const maxExtent = Math.max(1, ...points.map((point) => Math.max(Math.abs(point.re), Math.abs(point.im)))) * 1.15;
  const toX = (value) => width / 2 + (value / maxExtent) * (width / 2 - pad);
  const toY = (value) => height / 2 - (value / maxExtent) * (height / 2 - pad);

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fffefb";
  ctx.fillRect(0, 0, width, height);

  ctx.strokeStyle = "rgba(31, 26, 22, 0.12)";
  for (let i = -2; i <= 2; i += 1) {
    const grid = (i / 2) * maxExtent;
    ctx.beginPath();
    ctx.moveTo(toX(-maxExtent), toY(grid));
    ctx.lineTo(toX(maxExtent), toY(grid));
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(toX(grid), toY(-maxExtent));
    ctx.lineTo(toX(grid), toY(maxExtent));
    ctx.stroke();
  }

  ctx.strokeStyle = "#34291d";
  ctx.beginPath();
  ctx.moveTo(toX(-maxExtent), toY(0));
  ctx.lineTo(toX(maxExtent), toY(0));
  ctx.moveTo(toX(0), toY(-maxExtent));
  ctx.lineTo(toX(0), toY(maxExtent));
  ctx.stroke();

  ctx.fillStyle = "#66584b";
  ctx.font = "12px Segoe UI";
  ctx.fillText("In-phase", width - 84, height / 2 - 8);
  ctx.fillText("Quadrature", width / 2 + 8, 20);

  ctx.fillStyle = "rgba(15, 118, 110, 0.95)";
  references.forEach((point) => {
    ctx.beginPath();
    ctx.arc(toX(point.re), toY(point.im), 6, 0, Math.PI * 2);
    ctx.fill();
  });

  ctx.fillStyle = "rgba(180, 83, 9, 0.45)";
  received.forEach((point) => {
    ctx.beginPath();
    ctx.arc(toX(point.re), toY(point.im), 4, 0, Math.PI * 2);
    ctx.fill();
  });
}

function renderFixedConfig(config) {
  dom.fixedConfig.innerHTML = [
    `Scenario: ${config.scenario === "assignment2" ? "Assignment 2" : "Assignment 3"}`,
    `FFT Size: ${config.fftSize}`,
    `GI Ratio: ${config.giRatio}`,
    `Modulation: ${config.modulationOrder === 2 ? "QPSK" : "16QAM"}`,
    `Multipath: ${config.multiPath}`,
    `SNR: ${config.snrStart} to ${config.snrEnd} dB, step ${config.snrStep}`,
  ].join("<br>");
}

function renderLegend(series) {
  dom.legend.innerHTML = series.map((serie) => `
    <div class="legend-item">
      <span class="legend-swatch" style="background:${serie.color}"></span>
      <span>${serie.label}</span>
    </div>
  `).join("");
}

function renderLiveCards(series, iteration, targetIterations) {
  dom.liveCards.innerHTML = `
    <div class="summary-card">
      <div>Iteration Progress</div>
      <strong>${iteration} / ${targetIterations}</strong>
      <div>Curves are updated as batches finish.</div>
    </div>
    ${series.map((serie) => `
      <div class="summary-card">
        <div>${serie.label}</div>
        <strong>${Math.min(...serie.values.filter((value) => value > 0)).toExponential(2)}</strong>
        <div>Current best BER</div>
      </div>
    `).join("")}
  `;
}

function renderTable(result) {
  dom.resultsHeader.textContent = result.series.map((serie) => serie.label).join(" / ");
  dom.resultsTable.innerHTML = result.snrValues.map((snr, index) => {
    const values = result.series.map((serie) => {
      const value = serie.values[index];
      return value <= 0 ? "0" : (value < 0.001 ? value.toExponential(2) : value.toFixed(4));
    }).join(" / ");
    return `<tr><td>${snr}</td><td>${values}</td></tr>`;
  }).join("");
}

function drawChart(result) {
  const canvas = dom.chartCanvas;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  const pad = { left: 72, right: 24, top: 22, bottom: 52 };
  const plotWidth = width - pad.left - pad.right;
  const plotHeight = height - pad.top - pad.bottom;
  const xMin = result.snrValues[0];
  const xMax = result.snrValues[result.snrValues.length - 1];
  const xFor = (value) => xMax === xMin ? pad.left + plotWidth / 2 : pad.left + ((value - xMin) / (xMax - xMin)) * plotWidth;
  const yFor = (value) => {
    const clamped = Math.max(1e-4, Math.min(1, value || 1));
    const logValue = Math.log10(clamped);
    return pad.top + ((0 - logValue) / 4) * plotHeight;
  };

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fffefb";
  ctx.fillRect(0, 0, width, height);

  ctx.strokeStyle = "rgba(31, 26, 22, 0.12)";
  for (let i = 0; i <= 4; i += 1) {
    const y = pad.top + (plotHeight / 4) * i;
    ctx.beginPath();
    ctx.moveTo(pad.left, y);
    ctx.lineTo(width - pad.right, y);
    ctx.stroke();
  }

  ctx.strokeStyle = "#34291d";
  ctx.beginPath();
  ctx.moveTo(pad.left, pad.top);
  ctx.lineTo(pad.left, height - pad.bottom);
  ctx.lineTo(width - pad.right, height - pad.bottom);
  ctx.stroke();

  ctx.fillStyle = "#66584b";
  ctx.font = "12px Segoe UI";
  [1, 1e-1, 1e-2, 1e-3, 1e-4].forEach((tick) => ctx.fillText(tick.toExponential(0), 16, yFor(tick) + 4));
  result.snrValues.forEach((snr) => ctx.fillText(String(snr), xFor(snr) - 8, height - 20));
  ctx.fillText("SNR (dB)", width / 2 - 18, height - 8);

  result.series.forEach((serie) => {
    ctx.strokeStyle = serie.color;
    ctx.fillStyle = serie.color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    serie.values.forEach((value, index) => {
      const x = xFor(result.snrValues[index]);
      const y = yFor(value);
      if (index === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();
    serie.values.forEach((value, index) => {
      ctx.beginPath();
      ctx.arc(xFor(result.snrValues[index]), yFor(value), 4, 0, Math.PI * 2);
      ctx.fill();
    });
  });
}

function createProgressSeries(config) {
  const snrValues = createSnrRange(config.snrStart, config.snrEnd, config.snrStep);
  if (config.scenario === "assignment2") {
    return {
      snrValues,
      series: [{ key: "siso", label: "SISO-OFDM", color: palette.siso, values: Array(snrValues.length).fill(1) }],
      errors: { siso: Array(snrValues.length).fill(0) },
    };
  }

  return {
    snrValues,
    series: [
      { key: "sc", label: "SC", color: palette.sc, values: Array(snrValues.length).fill(1) },
      { key: "mrc", label: "MRC", color: palette.mrc, values: Array(snrValues.length).fill(1) },
      { key: "sskmSc", label: "SSKM-SC", color: palette.sskmSc, values: Array(snrValues.length).fill(1) },
      { key: "sskmMrc", label: "SSKM-MRC", color: palette.sskmMrc, values: Array(snrValues.length).fill(1) },
    ],
    errors: {
      sc: Array(snrValues.length).fill(0),
      mrc: Array(snrValues.length).fill(0),
      sskmSc: Array(snrValues.length).fill(0),
      sskmMrc: Array(snrValues.length).fill(0),
    },
  };
}

function runAssignment2Iteration(config, rng, state, snrIndex) {
  const dataSize = config.fftSize * config.modulationOrder;
  const giSize = Math.floor(config.fftSize * config.giRatio);
  const snr = state.snrValues[snrIndex];
  const data = randomBits(dataSize, rng);
  const modulated = baseMod(data, config.modulationOrder);
  const tx = addCyclicPrefix(normalizeIfft(modulated), giSize);
  const channel = rayleighChannel(config.multiPath, rng);
  const rx = addAwgn(convolve(tx, channel), snr, rng);
  const equalized = equalize(normalizeFft(removeCyclicPrefix(rx, giSize, config.fftSize)), channelResponse(channel, config.fftSize));
  state.errors.siso[snrIndex] += countBitErrors(data, baseDemod(equalized, config.modulationOrder).slice(0, data.length));
}

function runAssignment3Iteration(config, rng, state, snrIndex) {
  const giSize = Math.floor(config.fftSize * config.giRatio);
  const fftSizeSskm = config.fftSize * 2;
  const giSizeSskm = Math.floor(fftSizeSskm * config.giRatio);
  const dataSize = config.fftSize * config.modulationOrder;
  const dataSizeSskm = dataSize * 2;
  const snr = state.snrValues[snrIndex];

  const data = randomBits(dataSize, rng);
  const tx = addCyclicPrefix(normalizeIfft(baseMod(data, config.modulationOrder)), giSize);
  const h1 = rayleighChannel(config.multiPath, rng);
  const h2 = rayleighChannel(config.multiPath, rng);
  const rx1 = normalizeFft(removeCyclicPrefix(addAwgn(convolve(tx, h1), snr, rng), giSize, config.fftSize));
  const rx2 = normalizeFft(removeCyclicPrefix(addAwgn(convolve(tx, h2), snr, rng), giSize, config.fftSize));
  const H1 = channelResponse(h1, config.fftSize);
  const H2 = channelResponse(h2, config.fftSize);

  const sc = selectCombine(rx1, rx2, H1, H2);
  const mrc = mrcCombine(rx1, rx2, H1, H2);
  state.errors.sc[snrIndex] += countBitErrors(data, baseDemod(equalize(sc.y, sc.h), config.modulationOrder).slice(0, data.length));
  state.errors.mrc[snrIndex] += countBitErrors(data, baseDemod(equalize(mrc.y, mrc.h), config.modulationOrder).slice(0, data.length));

  const dataSskm = encodeSskm(data);
  const txSskm = addCyclicPrefix(normalizeIfft(baseMod(dataSskm, config.modulationOrder)), giSizeSskm);
  const rx1Sskm = normalizeFft(removeCyclicPrefix(addAwgn(convolve(txSskm, h1), snr, rng), giSizeSskm, fftSizeSskm));
  const rx2Sskm = normalizeFft(removeCyclicPrefix(addAwgn(convolve(txSskm, h2), snr, rng), giSizeSskm, fftSizeSskm));
  const H1Sskm = channelResponse(h1, fftSizeSskm);
  const H2Sskm = channelResponse(h2, fftSizeSskm);
  const scSskm = selectCombine(rx1Sskm, rx2Sskm, H1Sskm, H2Sskm);
  const mrcSskm = mrcCombine(rx1Sskm, rx2Sskm, H1Sskm, H2Sskm);

  state.errors.sskmSc[snrIndex] += countBitErrors(data, decodeSskm(baseDemod(equalize(scSskm.y, scSskm.h), config.modulationOrder).slice(0, dataSizeSskm)));
  state.errors.sskmMrc[snrIndex] += countBitErrors(data, decodeSskm(baseDemod(equalize(mrcSskm.y, mrcSskm.h), config.modulationOrder).slice(0, dataSizeSskm)));
}

function updateProgressValues(config, state, iteration) {
  const dataSize = config.fftSize * config.modulationOrder;
  state.series.forEach((serie) => {
    serie.values = state.errors[serie.key].map((error) => error / (dataSize * iteration));
  });
}

async function runBerProgress() {
  const token = ++runToken;
  const config = getBerConfig();
  const rng = new RNG(config.seed);
  const state = createProgressSeries(config);
  renderFixedConfig(config);
  renderLegend(state.series);

  const batchSize = Math.max(5, Math.floor(config.targetIterations / 25));
  for (let iteration = 1; iteration <= config.targetIterations; iteration += 1) {
    for (let snrIndex = 0; snrIndex < state.snrValues.length; snrIndex += 1) {
      if (config.scenario === "assignment2") runAssignment2Iteration(config, rng, state, snrIndex);
      else runAssignment3Iteration(config, rng, state, snrIndex);
    }

    if (iteration % batchSize === 0 || iteration === config.targetIterations) {
      if (token !== runToken) return;
      updateProgressValues(config, state, iteration);
      drawChart(state);
      renderLiveCards(state.series, iteration, config.targetIterations);
      renderTable(state);
      dom.statusText.textContent = `Running... ${iteration} / ${config.targetIterations}`;
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  }

  dom.statusText.textContent = `Done · ${config.targetIterations} iterations`;
}

function updateMode() {
  const isAssignment1 = dom.scenario.value === "assignment1";
  dom.assignment1Section.hidden = !isAssignment1;
  dom.berSection.hidden = isAssignment1;
  dom.berTopControls.hidden = isAssignment1;
  dom.presetButton.textContent = isAssignment1 ? "Load A1 Default" : "Load BER Default";
}

function applyDefaults() {
  if (dom.scenario.value === "assignment1") {
    dom.a1SampleCount.value = "256";
    dom.a1ChannelGain.value = "1";
    dom.a1ChannelPhase.value = "35";
    dom.a1NoiseSigma.value = "0.35";
    dom.seed.value = "20260318";
    updateAssignment1Labels();
  } else {
    dom.berTargetIterations.value = "500";
    updateBerTargetLabel();
    dom.seed.value = "20260318";
  }
}

async function runCurrentMode() {
  dom.runButton.disabled = true;
  updateMode();
  try {
    if (dom.scenario.value === "assignment1") {
      simulateAssignment1();
    } else {
      await runBerProgress();
    }
  } catch (error) {
    dom.statusText.textContent = error.message;
    dom.assignment1Status.textContent = error.message;
  } finally {
    dom.runButton.disabled = false;
  }
}

dom.runButton.addEventListener("click", runCurrentMode);
dom.presetButton.addEventListener("click", () => {
  applyDefaults();
  runCurrentMode();
});
dom.scenario.addEventListener("change", () => {
  updateMode();
  applyDefaults();
  runCurrentMode();
});
dom.berTargetIterations.addEventListener("input", updateBerTargetLabel);
[dom.a1SampleCount, dom.a1ChannelGain, dom.a1ChannelPhase, dom.a1NoiseSigma].forEach((input) => {
  input.addEventListener("input", () => {
    updateAssignment1Labels();
    if (dom.scenario.value === "assignment1") simulateAssignment1();
  });
});

updateAssignment1Labels();
updateBerTargetLabel();
updateMode();
runCurrentMode();
