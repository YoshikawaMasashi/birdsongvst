use crate::pi::{
    event::{Event, Triggered},
    types::AudioProcessor,
    utils::frequency_from_note_number,
};

fn bytes_to_f64_array(buf: Vec<u8>) -> Vec<f64> {
    let float_buf: Vec<f64> = buf
        .chunks(8)
        .map(|chunk| {
            let mut bytes_array: [u8; 8] = [0; 8];
            bytes_array.copy_from_slice(chunk);
            f64::from_le_bytes(bytes_array)
        })
        .collect();
    float_buf
}


fn ode(
    t: usize,
    v: [f64; 6],
    pi: &mut Vec<f64>,
    pb: &mut Vec<f64>,
    alpha: &Vec<f64>,
    beta: &Vec<f64>,
    gamma: f64,
    envelope: &Vec<f64>,
    c: f64,
    L: f64,
    r: f64,
    Ch: f64,
    MG: f64,
    MB: f64,
    RB: f64,
    Rh: f64,
    dt: f64,
    ovfs: usize,
) -> [f64; 6] {
    let [x, y, pout, i1, i2, i3] = v;
    let dv0 = y;
    let dv1 = (-alpha[t / ovfs] - beta[t / ovfs] * x - x * x * x + x * x) * gamma * gamma
        - (x * x * y + x * y) * gamma;

    let pbold = pb[t];
    let delayed_t = if t > ((L / c / dt) as usize) {t - ((L / c / dt) as usize)} else {0};
    pi[t] = (0.5 * envelope[t / ovfs]) * dv1 + pb[delayed_t];
    pb[t] = -r * pi[delayed_t];
    let pout = (1.0 - r) * pi[delayed_t];

    let dv2 = (pb[t] - pbold) / dt;
    let dv3 = i2;
    let dv4 = -(1.0 / Ch / MG) * i1 - Rh * (1.0 / MB + 1.0 / MG) * i2
        + (1.0 / MG / Ch + Rh * RB / MG / MB) * i3
        + (1.0 / MG) * dv2
        + (Rh * RB / MG / MB) * pout;
    let dv5 = -(MG / MB) * i2 - (Rh / MB) * i3 + (1.0 / MB) * pout;
    [dv0, dv1, dv2, dv3, dv4, dv5]
}

pub struct SineOscillator {
    phase: f64, // オシレータの進行具合。sin(x)のxを保持してる。
    freq: f64,  // 鳴らしている音の周波数。
    pitch: f64, // ピッチベンド幅を保持しておくフィールド。

    t: usize,

    v: [f64; 6],

    alpha: Vec<f64>,
    beta: Vec<f64>,
    envelope: Vec<f64>,
}

impl SineOscillator {
    pub fn new() -> Self {
        let alpha_bytes = std::include_bytes!("../../../rust_synth/data/alpha.raw");
        let alpha = bytes_to_f64_array(alpha_bytes.to_vec());
        let beta_bytes = std::include_bytes!("../../../rust_synth/data/beta.raw");
        let beta = bytes_to_f64_array(beta_bytes.to_vec());
        let envelope_bytes = std::include_bytes!("../../../rust_synth/data/envelope.raw");
        let envelope = bytes_to_f64_array(envelope_bytes.to_vec());
        Self {
            phase: 0.0,
            freq: 440.0,
            pitch: 0.0,

            t: 0,
            v: [1e-4 * 1e2, 1e-4 * 1e1, 1e-4, 1e-4, 1e-4, 1e-4],

            alpha,
            beta,
            envelope,
        }
    }
}

impl Triggered for SineOscillator {
    fn trigger(&mut self, event: &Event) {
        match event {
            Event::NoteOn { note, velocity: _ } => {
                self.t = 0;
            }
            Event::NoteOff { note: _ } => {}
            Event::PitchBend { ratio } => {}
        }
    }
}

impl AudioProcessor<f64> for SineOscillator {
    fn process(&mut self, sample_rate: f64) -> f64 {
        0.0
    }
}
