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

pub struct SineOscillator {
    phase: f64, // オシレータの進行具合。sin(x)のxを保持してる。
    freq: f64,  // 鳴らしている音の周波数。
    pitch: f64, // ピッチベンド幅を保持しておくフィールド。

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
                self.freq = frequency_from_note_number(*note);
            }
            Event::NoteOff { note: _ } => {}
            Event::PitchBend { ratio } => {
                self.pitch = *ratio;
            }
        }
    }
}

impl AudioProcessor<f64> for SineOscillator {
    fn process(&mut self, sample_rate: f64) -> f64 {
        let phase_diff = (self.freq * self.pitch) * 2.0 * std::f64::consts::PI / sample_rate;
        self.phase += phase_diff;

        self.phase.sin()
    }
}
