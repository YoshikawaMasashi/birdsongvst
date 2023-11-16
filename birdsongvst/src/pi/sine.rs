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
    fs: usize,
    gamma: f64,
    ovfs: usize,
    
    t: usize,

    v: [f64; 6],

    alpha: Vec<f64>,
    beta: Vec<f64>,
    envelope: Vec<f64>,
    pi: Vec<f64>,
    pb: Vec<f64>,
}

impl SineOscillator {
    pub fn new() -> Self {
        let alpha_bytes = std::include_bytes!("../../../rust_synth/data/alpha.raw");
        let alpha = bytes_to_f64_array(alpha_bytes.to_vec());
        let beta_bytes = std::include_bytes!("../../../rust_synth/data/beta.raw");
        let beta = bytes_to_f64_array(beta_bytes.to_vec());
        let envelope_bytes = std::include_bytes!("../../../rust_synth/data/envelope.raw");
        let envelope = bytes_to_f64_array(envelope_bytes.to_vec());
        let ovfs = 20;
        let tmax = alpha.len() * ovfs - 1;
        Self {
            fs: 44100,
            gamma: 40000.0,
            ovfs,
            t: 0,
            v: [1e-4 * 1e2, 1e-4 * 1e1, 1e-4, 1e-4, 1e-4, 1e-4],

            alpha,
            beta,
            envelope,
            pi: vec![0.0; tmax],
            pb: vec![0.0; tmax],
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
        let (c, L, r, Ch, MG, MB, RB, Rh) = (
            343.0, 0.025, 0.65, 1.43e-10, 20.0, 10000.0, 5000000.0, 24000.0,
        );
        let tmax = self.alpha.len() * self.ovfs - 1;
        let dt = 1.0 / (self.ovfs * self.fs) as f64;

        for i in 0..self.ovfs {
            if self.t >= tmax {
                return 0.0;
            }

            let k1 = ode(
                self.t, self.v, &mut self.pi, &mut self.pb, &self.alpha, &self.beta, self.gamma, &self.envelope, c, L, r, Ch, MG, MB, RB, Rh,
                dt, self.ovfs,
            );
    
            let v2_vec: Vec<f64> = self.v
                .iter()
                .zip(k1.iter())
                .map(|(v_, k1_)| v_ + dt / 2.0 * k1_)
                .collect();
            let mut v2: [f64; 6] = [0.0; 6];
            v2.copy_from_slice(&v2_vec);
            let k2 = ode(
                self.t, v2, &mut self.pi, &mut self.pb, &self.alpha, &self.beta, self.gamma, &self.envelope, c, L, r, Ch, MG, MB, RB, Rh,
                dt, self.ovfs,
            );
    
            let v3_vec: Vec<f64> = self.v
                .iter()
                .zip(k2.iter())
                .map(|(v_, k2_)| v_ + dt / 2.0 * k2_)
                .collect();
            let mut v3: [f64; 6] = [0.0; 6];
            v3.copy_from_slice(&v3_vec);
            let k3 = ode(
                self.t, v3, &mut self.pi, &mut self.pb, &self.alpha, &self.beta, self.gamma, &self.envelope, c, L, r, Ch, MG, MB, RB, Rh,
                dt, self.ovfs,
            );
    
            let v4_vec: Vec<f64> = self.v
                .iter()
                .zip(k3.iter())
                .map(|(v_, k3_)| v_ + dt * k3_)
                .collect();
            let mut v4: [f64; 6] = [0.0; 6];
            v4.copy_from_slice(&v4_vec);
            let k4 = ode(
                self.t, v4, &mut self.pi, &mut self.pb, &self.alpha, &self.beta, self.gamma, &self.envelope, c, L, r, Ch, MG, MB, RB, Rh,
                dt, self.ovfs,
            );
    
            self.v = [
                self.v[0] + dt * (2.0 * (k2[0] + k3[0]) + k1[0] + k4[0]) / 6.0,
                self.v[1] + dt * (2.0 * (k2[1] + k3[1]) + k1[1] + k4[1]) / 6.0,
                self.v[2] + dt * (2.0 * (k2[2] + k3[2]) + k1[2] + k4[2]) / 6.0,
                self.v[3] + dt * (2.0 * (k2[3] + k3[3]) + k1[3] + k4[3]) / 6.0,
                self.v[4] + dt * (2.0 * (k2[4] + k3[4]) + k1[4] + k4[4]) / 6.0,
                self.v[5] + dt * (2.0 * (k2[5] + k3[5]) + k1[5] + k4[5]) / 6.0,
            ];

            self.t += 1;
        }
        // TODO: ノーマライズ
        RB * self.v[5]
    }
}
