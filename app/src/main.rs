//! Desktop front end: history tape, keypad, graphing, function catalog.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use cas::numeric::{self, evalf};
use cas::demos::{DEMOS, StepResult};
use cas::{Angle, CATALOG, Calc, Expr, Figure, Notation, Shape};
use eframe::egui::{self, Align, Color32, Key, Layout, RichText};
use egui_plot::{Arrows, Legend, Line, LineStyle, MarkerShape, Plot, PlotBounds, PlotPoint, PlotPoints, Points, Polygon, Text};

#[derive(PartialEq, Clone, Copy)]
enum Tab {
    Calc,
    Graph,
    Demos,
    Catalog,
}

#[derive(PartialEq, Clone, Copy)]
enum Kind {
    Function,
    Parametric,
    Polar,
}

struct Entry {
    input: String,
    text: String,
    approx: Option<String>,
    note: Option<String>,
    error: bool,
    figure: Option<Figure>,
    show_figure: bool,
}

struct PlotFn {
    kind: Kind,
    src: String,
    src2: String,
    on: bool,
    /// (source it was compiled from, result)
    compiled: Option<(String, Result<(Expr, Option<Expr>), String>)>,
}

impl PlotFn {
    fn new(src: &str) -> PlotFn {
        PlotFn { kind: Kind::Function, src: src.into(), src2: String::new(), on: true, compiled: None }
    }
}

struct App {
    calc: Calc,
    input: String,
    history: Vec<Entry>,
    recall: Option<usize>,
    tab: Tab,
    second: bool,
    refocus: bool,
    fns: Vec<PlotFn>,
    t_range: (f64, f64),
    bounds: PlotBounds,
    marks: Vec<([f64; 2], String)>,
    equal_axes: bool,
    reset_view: bool,
    filter: String,
    demo: usize,
    demo_steps: Vec<StepResult>,
}

const COLORS: [Color32; 6] = [
    Color32::from_rgb(66, 135, 245),
    Color32::from_rgb(235, 87, 87),
    Color32::from_rgb(39, 174, 96),
    Color32::from_rgb(242, 153, 74),
    Color32::from_rgb(155, 81, 224),
    Color32::from_rgb(45, 185, 195),
];

/// (label, inserted text, 2nd label, 2nd inserted text)
type KeyDef = (&'static str, &'static str, &'static str, &'static str);
const KEYS: &[&[KeyDef]] = &[
    &[("solve", "solve(", "csolve", "csolve("), ("d/dx", "diff(", "∫", "integrate("), ("lim", "limit(", "Σ", "sum("), ("factor", "factor(", "expand", "expand("), ("simp", "simplify(", "taylor", "taylor("), ("x", "x", "y", "y"), ("=", "=", ":=", ":=")],
    &[("sin", "sin(", "asin", "asin("), ("cos", "cos(", "acos", "acos("), ("tan", "tan(", "atan", "atan("), ("ln", "ln(", "e^", "e^("), ("log", "log(", "10^", "10^("), ("√", "sqrt(", "∛", "cbrt("), ("^", "^", "x²", "^2")],
    &[("π", "π", "∞", "∞"), ("e", "e", "i", "i"), ("(", "(", "[", "["), (")", ")", "]", "]"), (",", ",", ";", ";"), ("abs", "abs(", "n!", "!"), ("÷", "/", "1/x", "^-1")],
    &[("ans", "ans", "det", "det("), ("{", "{", "<", "<"), ("7", "7", "sinh", "sinh("), ("8", "8", "cosh", "cosh("), ("9", "9", "tanh", "tanh("), ("}", "}", ">", ">"), ("×", "*", "nCr", "ncr(")],
    &[("EE", "E", "inv", "inv("), ("|", " | ", "mean", "mean("), ("4", "4", "", ""), ("5", "5", "", ""), ("6", "6", "", ""), ("t", "t", "n", "n"), ("−", "-", "nPr", "npr(")],
    &[("(−)", "-", "mod", "mod("), ("0", "0", "", ""), ("1", "1", "", ""), ("2", "2", "", ""), ("3", "3", "", ""), (".", ".", "", ""), ("+", "+", "gcd", "gcd(")],
];

impl App {
    fn new() -> App {
        App {
            calc: Calc::new(),
            input: String::new(),
            history: Vec::new(),
            recall: None,
            tab: Tab::Calc,
            second: false,
            refocus: true,
            fns: vec![PlotFn::new("sin(x)"), PlotFn::new("")],
            t_range: (0.0, std::f64::consts::TAU),
            bounds: PlotBounds::from_min_max([-10.0, -10.0], [10.0, 10.0]),
            marks: Vec::new(),
            equal_axes: false,
            reset_view: false,
            filter: String::new(),
            demo: 0,
            demo_steps: cas::demos::run(&DEMOS[0]),
        }
    }

    fn submit(&mut self) {
        let src = self.input.trim().to_string();
        if src.is_empty() {
            return;
        }
        let s = &mut self.calc.settings;
        let mode = match src.as_str() {
            "deg" => Some(s.angle = Angle::Deg),
            "rad" => Some(s.angle = Angle::Rad),
            "grad" => Some(s.angle = Angle::Grad),
            "exact" => Some(s.approx = false),
            "approx" => Some(s.approx = true),
            "real" => Some(s.complex = false),
            "complex" => Some(s.complex = true),
            _ => None,
        };
        if src == "clear" {
            self.history.clear();
        } else if mode.is_some() {
            self.history.push(Entry { input: src.clone(), text: format!("mode: {src}"), approx: None, note: None, error: false, figure: None, show_figure: false });
        } else {
            let e = match self.calc.eval_line(&src) {
                Ok(o) => Entry { input: src, text: o.text, approx: o.approx, note: o.note, error: false, show_figure: o.figure.is_some(), figure: o.figure },
                Err(e) => Entry { input: src, text: e, approx: None, note: None, error: true, figure: None, show_figure: false },
            };
            // only the newest figure stays open, so the tape stays short
            for old in &mut self.history {
                old.show_figure = false;
            }
            self.history.push(e);
        }
        self.input.clear();
        self.recall = None;
        self.refocus = true;
        // stored variables may change what the graph functions mean
        for f in &mut self.fns {
            f.compiled = None;
        }
    }

    fn mode_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut self.tab, Tab::Calc, "Calculator");
            ui.selectable_value(&mut self.tab, Tab::Graph, "Graph");
            ui.selectable_value(&mut self.tab, Tab::Demos, "Demos");
            ui.selectable_value(&mut self.tab, Tab::Catalog, "Catalog");
            ui.separator();
            let s = &mut self.calc.settings;
            let before = s.angle;
            ui.selectable_value(&mut s.angle, Angle::Rad, "RAD");
            ui.selectable_value(&mut s.angle, Angle::Deg, "DEG");
            ui.selectable_value(&mut s.angle, Angle::Grad, "GRAD");
            if before != s.angle {
                for f in &mut self.fns {
                    f.compiled = None;
                }
            }
            ui.separator();
            ui.selectable_value(&mut s.approx, false, "Exact");
            ui.selectable_value(&mut s.approx, true, "Approx");
            ui.separator();
            ui.selectable_value(&mut s.complex, false, "Real");
            ui.selectable_value(&mut s.complex, true, "Complex");
            ui.separator();
            ui.selectable_value(&mut s.fmt.notation, Notation::Normal, "Norm");
            ui.selectable_value(&mut s.fmt.notation, Notation::Sci, "Sci");
            ui.selectable_value(&mut s.fmt.notation, Notation::Eng, "Eng");
            ui.separator();
            let mut fixed = s.fmt.fix.is_some();
            if ui.checkbox(&mut fixed, "Fix").changed() {
                s.fmt.fix = fixed.then_some(2);
            }
            match &mut s.fmt.fix {
                Some(n) => {
                    ui.add(egui::DragValue::new(n).range(0..=12));
                }
                None => {
                    ui.label("digits");
                    ui.add(egui::DragValue::new(&mut s.fmt.digits).range(3..=15));
                }
            }
        });
    }

    fn keypad(&mut self, ui: &mut egui::Ui) {
        let cols = KEYS[0].len() as f32 + 1.0;
        let w = ((ui.available_width() - 4.0 * cols) / cols).clamp(34.0, 150.0);
        let size = [w, 34.0];
        let mut insert: Option<&str> = None;
        let mut enter = false;
        egui::Grid::new("keys").spacing([4.0, 4.0]).show(ui, |ui| {
            for (r, row) in KEYS.iter().enumerate() {
                // left utility column
                match r {
                    0 => {
                        let b = egui::Button::new(RichText::new("2nd").strong()).selected(self.second);
                        if ui.add_sized(size, b).clicked() {
                            self.second = !self.second;
                        }
                    }
                    1 => {
                        if ui.add_sized(size, egui::Button::new("DEL")).clicked() {
                            self.input.pop();
                            self.refocus = true;
                        }
                    }
                    2 => {
                        if ui.add_sized(size, egui::Button::new("CLR")).clicked() {
                            self.input.clear();
                            self.refocus = true;
                        }
                    }
                    5 => {
                        if ui.add_sized(size, egui::Button::new(RichText::new("ENTER").strong())).clicked() {
                            enter = true;
                        }
                    }
                    3 => {
                        if ui.add_sized(size, egui::Button::new("≈")).on_hover_text("decimal value of the last answer").clicked() {
                            insert = Some("approx(ans)");
                            enter = true;
                        }
                    }
                    _ => {
                        if ui.add_sized(size, egui::Button::new(RichText::new("k").size(15.0))).clicked() {
                            insert = Some("k");
                        }
                    }
                }
                for (label, text, label2, text2) in row.iter() {
                    let (l, t) = if self.second && !label2.is_empty() { (*label2, *text2) } else { (*label, *text) };
                    let mut b = egui::Button::new(RichText::new(l).size(15.0));
                    if l.chars().all(|c| c.is_ascii_digit() || c == '.') {
                        b = b.fill(ui.visuals().extreme_bg_color);
                    }
                    let resp = ui.add_sized(size, b);
                    let resp = if !self.second && !label2.is_empty() { resp.on_hover_text(format!("2nd: {label2}")) } else { resp };
                    if resp.clicked() {
                        insert = Some(t);
                    }
                }
                ui.end_row();
            }
        });
        if let Some(t) = insert {
            self.input.push_str(t);
            self.second = false;
            self.refocus = true;
        }
        if enter {
            self.submit();
        }
    }

    fn calc_tab(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("entry").show(ui, |ui| {
            ui.add_space(6.0);
            let edit = egui::TextEdit::singleline(&mut self.input).desired_width(f32::INFINITY).font(egui::FontId::monospace(18.0)).hint_text("type an expression, Enter to evaluate");
            let resp = ui.add(edit);
            if resp.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                self.submit();
            }
            if resp.has_focus() && !self.history.is_empty() {
                let (up, down) = ui.input(|i| (i.key_pressed(Key::ArrowUp), i.key_pressed(Key::ArrowDown)));
                if up || down {
                    let n = self.history.len();
                    let pos = match (self.recall, up) {
                        (None, true) => Some(n - 1),
                        (Some(p), true) => Some(p.saturating_sub(1)),
                        (Some(p), false) if p + 1 < n => Some(p + 1),
                        _ => None,
                    };
                    self.recall = pos;
                    self.input = pos.map(|p| self.history[p].input.clone()).unwrap_or_default();
                }
            }
            if self.refocus {
                resp.request_focus();
                self.refocus = false;
            }
            ui.horizontal(|ui| {
                ui.label(RichText::new("show last answer as").small().weak());
                for (label, cmd) in [("rad", "torad(ans)"), ("DD", "todd(ans)"), ("DMS", "todms(ans)")] {
                    if ui.small_button(label).clicked() {
                        self.input = cmd.to_string();
                        self.submit();
                    }
                }
            });
            ui.add_space(4.0);
            self.keypad(ui);
            ui.add_space(4.0);
        });
        egui::CentralPanel::default().show(ui, |ui| {
            let mut paste: Option<String> = None;
            egui::ScrollArea::vertical().auto_shrink([false, false]).stick_to_bottom(true).show(ui, |ui| {
                if self.history.is_empty() {
                    ui.add_space(12.0);
                    ui.label(RichText::new("Examples (click to try)").weak());
                    for ex in ["solve(x^2 - 5x + 6 = 0, x)", "integrate(x*e^x, x)", "diff(sin(x)^2, x)", "limit(sin(x)/x, x, 0)", "factor(x^4 - 16)", "det([[1, 2], [3, 4]])", "sum(1/k^2, k, 1, inf)", "convert(60, mph, kph)"] {
                        if ui.link(RichText::new(ex).monospace()).clicked() {
                            paste = Some(ex.to_string());
                        }
                    }
                }
                let mut to_graph: Option<String> = None;
                for (idx, e) in self.history.iter_mut().enumerate() {
                    ui.add_space(4.0);
                    if ui.add(egui::Label::new(RichText::new(&e.input).monospace().size(15.0).weak()).sense(egui::Sense::click())).on_hover_text("click to edit again").clicked() {
                        paste = Some(e.input.clone());
                    }
                    ui.with_layout(Layout::top_down(Align::Max), |ui| {
                        let color = if e.error { ui.visuals().error_fg_color } else { ui.visuals().strong_text_color() };
                        let main = egui::Label::new(RichText::new(&e.text).monospace().size(19.0).color(color)).wrap().sense(egui::Sense::click());
                        if ui.add(main).on_hover_text("click to insert").clicked() && !e.error {
                            paste = Some(e.text.clone());
                        }
                        if let Some(a) = &e.approx {
                            let lbl = egui::Label::new(RichText::new(format!("≈ {a}")).monospace().size(15.0)).wrap().sense(egui::Sense::click());
                            if ui.add(lbl).clicked() {
                                paste = Some(a.clone());
                            }
                        }
                        if let Some(n) = &e.note {
                            ui.label(RichText::new(n).small().weak());
                        }
                    });
                    if let Some(fig) = &e.figure {
                        ui.horizontal(|ui| {
                            if ui.small_button(if e.show_figure { "hide figure" } else { "show figure" }).clicked() {
                                e.show_figure = !e.show_figure;
                            }
                            if let Some(src) = &fig.graph_src {
                                if ui.small_button("open in Graph").clicked() {
                                    to_graph = Some(src.clone());
                                }
                            }
                        });
                        if e.show_figure {
                            draw_figure(ui, ("fig", idx), fig, 280.0);
                        }
                    }
                    ui.separator();
                }
                if let Some(src) = to_graph {
                    match self.fns.iter_mut().find(|f| f.src.trim().is_empty()) {
                        Some(slot) => *slot = PlotFn::new(&src),
                        None => self.fns.push(PlotFn::new(&src)),
                    }
                    self.tab = Tab::Graph;
                }
            });
            if let Some(p) = paste {
                self.input.push_str(&p);
                self.refocus = true;
            }
        });
    }

    fn demos_tab(&mut self, ui: &mut egui::Ui) {
        egui::Panel::left("demo_list").min_size(230.0).show(ui, |ui| {
            egui::ScrollArea::vertical().id_salt("demo_scroll").show(ui, |ui| {
                let mut group = "";
                for (i, d) in DEMOS.iter().enumerate() {
                    if d.group != group {
                        group = d.group;
                        ui.add_space(6.0);
                        ui.label(RichText::new(group).strong());
                    }
                    if ui.selectable_label(self.demo == i, d.title).clicked() && self.demo != i {
                        self.demo = i;
                        self.demo_steps = cas::demos::run(d);
                    }
                }
            });
        });
        egui::CentralPanel::default().show(ui, |ui| {
            let d = &DEMOS[self.demo];
            let mut try_it: Option<String> = None;
            let mut run_all = false;
            egui::ScrollArea::vertical().auto_shrink([false, false]).id_salt(("demo_body", self.demo)).show(ui, |ui| {
                ui.heading(d.title);
                ui.label(d.about);
                ui.add_space(4.0);
                if ui.button("Run every step in the calculator").clicked() {
                    run_all = true;
                }
                ui.separator();
                for (i, step) in self.demo_steps.iter().enumerate() {
                    if !step.note.is_empty() {
                        ui.add_space(6.0);
                        ui.label(step.note);
                    }
                    ui.horizontal_wrapped(|ui| {
                        ui.label(RichText::new("type:").small().weak());
                        if ui.link(RichText::new(step.input).monospace().size(15.0)).on_hover_text("click to put this in the calculator").clicked() {
                            try_it = Some(step.input.to_string());
                        }
                    });
                    match &step.output {
                        Ok(o) => {
                            ui.add(egui::Label::new(RichText::new(o).monospace().size(16.0).color(ui.visuals().strong_text_color())).wrap());
                        }
                        Err(e) => {
                            ui.colored_label(ui.visuals().error_fg_color, e);
                        }
                    }
                    if let Some(fig) = &step.figure {
                        draw_figure(ui, ("demo_fig", self.demo, i), fig, 300.0);
                    }
                }
            });
            if run_all {
                for (_, input) in d.steps {
                    match *input {
                        "deg" => self.calc.settings.angle = Angle::Deg,
                        "rad" => self.calc.settings.angle = Angle::Rad,
                        src => {
                            self.input = src.to_string();
                            self.submit();
                        }
                    }
                }
                self.tab = Tab::Calc;
            }
            if let Some(t) = try_it {
                self.input = t;
                self.tab = Tab::Calc;
                self.refocus = true;
            }
        });
    }

    fn compile(&mut self) {
        for f in &mut self.fns {
            let key = format!("{}\u{1}{}", f.src, f.src2);
            if f.compiled.as_ref().is_some_and(|(k, _)| *k == key) {
                continue;
            }
            let var = if f.kind == Kind::Function { "x" } else { "t" };
            let r = if f.src.trim().is_empty() {
                Err(String::new())
            } else {
                self.calc.plot_fn(&f.src, var).and_then(|a| {
                    if f.kind == Kind::Parametric {
                        if f.src2.trim().is_empty() {
                            return Err("enter y(t)".into());
                        }
                        Ok((a, Some(self.calc.plot_fn(&f.src2, var)?)))
                    } else {
                        Ok((a, None))
                    }
                })
            };
            f.compiled = Some((key, r));
        }
    }

    /// Curves as lists of unbroken segments.
    fn sample(&self, f: &PlotFn) -> Vec<Vec<[f64; 2]>> {
        let Some((_, Ok((a, b)))) = &f.compiled else { return vec![] };
        let (min, max) = (self.bounds.min(), self.bounds.max());
        let n = 1600;
        let pts: Vec<[f64; 2]> = match f.kind {
            Kind::Function => (0..=n)
                .map(|i| {
                    let x = min[0] + (max[0] - min[0]) * i as f64 / n as f64;
                    [x, evalf(a, "x", x)]
                })
                .collect(),
            Kind::Parametric | Kind::Polar => (0..=n)
                .map(|i| {
                    let t = self.t_range.0 + (self.t_range.1 - self.t_range.0) * i as f64 / n as f64;
                    match b {
                        Some(b) => [evalf(a, "t", t), evalf(b, "t", t)],
                        None => {
                            let r = evalf(a, "t", t);
                            [r * t.cos(), r * t.sin()]
                        }
                    }
                })
                .collect(),
        };
        // break at undefined points and at jumps across an asymptote
        let jump = 4.0 * (max[1] - min[1]).abs().max(1e-9);
        let mut out: Vec<Vec<[f64; 2]>> = vec![vec![]];
        for p in pts {
            let bad = !p[0].is_finite() || !p[1].is_finite();
            let far = out.last().unwrap().last().is_some_and(|q: &[f64; 2]| f.kind == Kind::Function && (q[1] - p[1]).abs() > jump);
            if bad || far {
                out.push(vec![]);
            }
            if !bad {
                out.last_mut().unwrap().push(p);
            }
        }
        out.retain(|s| s.len() > 1);
        out
    }

    fn analyze(&mut self, what: &str) {
        let (lo, hi) = (self.bounds.min()[0], self.bounds.max()[0]);
        let fs: Vec<Expr> = self.fns.iter().filter(|f| f.on && f.kind == Kind::Function).filter_map(|f| f.compiled.as_ref().and_then(|(_, r)| r.as_ref().ok().map(|(a, _)| a.clone()))).collect();
        let lbl = |x: f64, y: f64| format!("({}, {})", trim(x), trim(y));
        match what {
            "zeros" => {
                for f in &fs {
                    for x in numeric::find_roots(&|v| evalf(f, "x", v), lo, hi, 3000) {
                        self.marks.push(([x, 0.0], lbl(x, 0.0)));
                    }
                }
            }
            "extrema" => {
                for f in &fs {
                    let g = |v: f64| evalf(f, "x", v);
                    for x in numeric::find_roots(&|v| numeric::nderiv(&g, v), lo, hi, 3000) {
                        let h = 1e-3 * (1.0 + x.abs());
                        let curv = g(x + h) - 2.0 * g(x) + g(x - h);
                        if curv.abs() > 1e-12 && g(x).is_finite() {
                            self.marks.push(([x, g(x)], format!("{} {}", if curv > 0.0 { "min" } else { "max" }, lbl(x, g(x)))));
                        }
                    }
                }
            }
            _ => {
                for i in 0..fs.len() {
                    for j in i + 1..fs.len() {
                        for x in numeric::find_roots(&|v| evalf(&fs[i], "x", v) - evalf(&fs[j], "x", v), lo, hi, 3000) {
                            let y = evalf(&fs[i], "x", x);
                            self.marks.push(([x, y], lbl(x, y)));
                        }
                    }
                }
            }
        }
    }

    fn graph_tab(&mut self, ui: &mut egui::Ui) {
        self.compile();
        egui::Panel::left("fns").min_size(250.0).show(ui, |ui| {
            ui.add_space(4.0);
            let mut remove = None;
            for (i, f) in self.fns.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.checkbox(&mut f.on, "");
                    ui.colored_label(COLORS[i % COLORS.len()], "■");
                    let before = f.kind;
                    egui::ComboBox::from_id_salt(("kind", i))
                        .width(64.0)
                        .selected_text(match f.kind {
                            Kind::Function => "y(x)",
                            Kind::Parametric => "x(t), y(t)",
                            Kind::Polar => "r(t)",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut f.kind, Kind::Function, "y(x)");
                            ui.selectable_value(&mut f.kind, Kind::Parametric, "x(t), y(t)");
                            ui.selectable_value(&mut f.kind, Kind::Polar, "r(t)  polar");
                        });
                    if before != f.kind {
                        f.compiled = None;
                    }
                    if ui.small_button("✖").clicked() {
                        remove = Some(i);
                    }
                });
                let hint = match f.kind {
                    Kind::Function => "y = f(x)",
                    Kind::Parametric => "x(t)",
                    Kind::Polar => "r(t)",
                };
                ui.add(egui::TextEdit::singleline(&mut f.src).desired_width(f32::INFINITY).font(egui::FontId::monospace(15.0)).hint_text(hint));
                if f.kind == Kind::Parametric {
                    ui.add(egui::TextEdit::singleline(&mut f.src2).desired_width(f32::INFINITY).font(egui::FontId::monospace(15.0)).hint_text("y(t)"));
                }
                if let Some((_, Err(e))) = &f.compiled {
                    if !e.is_empty() {
                        ui.colored_label(ui.visuals().error_fg_color, e);
                    }
                }
                ui.add_space(4.0);
            }
            if let Some(i) = remove {
                self.fns.remove(i);
            }
            if ui.button("+ add function").clicked() {
                self.fns.push(PlotFn::new(""));
            }
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("t from");
                ui.add(egui::DragValue::new(&mut self.t_range.0).speed(0.1));
                ui.label("to");
                ui.add(egui::DragValue::new(&mut self.t_range.1).speed(0.1));
            });
            ui.checkbox(&mut self.equal_axes, "equal axis scale");
            ui.separator();
            ui.label("Analyze (visible window)");
            ui.horizontal_wrapped(|ui| {
                for (label, what) in [("Zeros", "zeros"), ("Min / Max", "extrema"), ("Intersections", "cross")] {
                    if ui.button(label).clicked() {
                        self.analyze(what);
                    }
                }
            });
            ui.horizontal(|ui| {
                if ui.button("Clear marks").clicked() {
                    self.marks.clear();
                }
                if ui.button("Reset view").clicked() {
                    self.reset_view = true;
                }
            });
            ui.add_space(6.0);
            ui.label(RichText::new("drag to pan, scroll to zoom, right-drag to box zoom, hover a curve to trace").small().weak());
        });
        egui::CentralPanel::default().show(ui, |ui| {
            let mut plot = Plot::new("plot").legend(Legend::default()).auto_bounds(egui::Vec2b::FALSE).default_x_bounds(-10.0, 10.0).default_y_bounds(-10.0, 10.0);
            if self.equal_axes {
                plot = plot.data_aspect(1.0);
            }
            if std::mem::take(&mut self.reset_view) {
                plot = plot.reset();
            }
            let curves: Vec<(usize, String, Vec<Vec<[f64; 2]>>)> = self.fns.iter().enumerate().filter(|(_, f)| f.on).map(|(i, f)| (i, if f.kind == Kind::Parametric { format!("({}, {})", f.src, f.src2) } else { f.src.clone() }, self.sample(f))).collect();
            let marks = self.marks.clone();
            let text_color = ui.visuals().strong_text_color();
            let resp = plot.show(ui, |pui| {
                for (i, name, segs) in curves {
                    for seg in segs {
                        pui.line(Line::new(name.clone(), PlotPoints::from(seg)).color(COLORS[i % COLORS.len()]).width(2.0));
                    }
                }
                if !marks.is_empty() {
                    pui.points(Points::new("marks", marks.iter().map(|m| m.0).collect::<Vec<_>>()).radius(4.5).shape(MarkerShape::Circle).color(text_color));
                    let dy = pui.plot_bounds().height() * 0.035;
                    for (p, label) in &marks {
                        pui.text(Text::new("mark", PlotPoint::new(p[0], p[1] + dy), label.clone()).color(text_color));
                    }
                }
                pui.plot_bounds()
            });
            self.bounds = resp.inner;
        });
    }

    fn catalog_tab(&mut self, ui: &mut egui::Ui) {
        egui::CentralPanel::default().show(ui, |ui| {
            ui.add_space(4.0);
            ui.add(egui::TextEdit::singleline(&mut self.filter).desired_width(f32::INFINITY).hint_text("search functions"));
            ui.add_space(4.0);
            let needle = self.filter.to_lowercase();
            let mut pick: Option<String> = None;
            let mut del: Option<String> = None;
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                let mut names: Vec<&String> = self.calc.vars.keys().chain(self.calc.funcs.keys()).collect();
                names.sort();
                if !names.is_empty() {
                    ui.label(RichText::new("Stored").strong());
                    for n in names {
                        ui.horizontal(|ui| {
                            let shown = match (self.calc.vars.get(n), self.calc.funcs.get(n)) {
                                (Some(v), _) => format!("{} = {}", n, self.calc.show(v)),
                                (_, Some((p, _))) => format!("{}({})", n, p.join(", ")),
                                _ => n.clone(),
                            };
                            if ui.small_button("✖").clicked() {
                                del = Some(n.clone());
                            }
                            if ui.link(RichText::new(shown).monospace()).clicked() {
                                pick = Some(n.clone());
                            }
                        });
                    }
                    ui.separator();
                }
                for (syntax, desc) in CATALOG.iter().filter(|(s, d)| needle.is_empty() || s.to_lowercase().contains(&needle) || d.to_lowercase().contains(&needle)) {
                    ui.horizontal_wrapped(|ui| {
                        if ui.link(RichText::new(*syntax).monospace().size(15.0)).clicked() {
                            // insert the call up to its opening parenthesis
                            pick = Some(match syntax.find('(') {
                                Some(k) if syntax.ends_with(')') && !syntax.starts_with('[') => syntax[..=k].to_string(),
                                _ => syntax.to_string(),
                            });
                        }
                        ui.label(RichText::new(*desc).weak());
                    });
                }
            });
            if let Some(n) = del {
                self.calc.vars.remove(&n);
                self.calc.funcs.remove(&n);
            }
            if let Some(p) = pick {
                self.input.push_str(&p);
                self.tab = Tab::Calc;
                self.refocus = true;
            }
        });
    }
}

fn figure_color(i: u8, ui: &egui::Ui) -> Color32 {
    if i >= 6 { ui.visuals().weak_text_color() } else { COLORS[i as usize % COLORS.len()] }
}

/// Draw a figure made of plain shapes.
fn draw_figure(ui: &mut egui::Ui, id: impl std::hash::Hash + std::fmt::Debug, fig: &Figure, height: f32) {
    ui.label(RichText::new(&fig.title).strong());
    let mut plot = Plot::new(id).height(height).legend(Legend::default()).allow_scroll(false).show_axes([!fig.equal, !fig.equal]).show_grid([!fig.equal, !fig.equal]);
    if fig.equal || fig.square {
        // a square frame keeps the whole figure in view
        let side = height.max(400.0).min(ui.available_width());
        plot = plot.width(side).height(side);
    }
    if fig.equal {
        plot = plot.data_aspect(1.0);
    }
    if let Some(b) = fig.bounds {
        plot = plot.auto_bounds(egui::Vec2b::FALSE).default_x_bounds(b[0], b[1]).default_y_bounds(b[2], b[3]);
    }
    let text_color = ui.visuals().strong_text_color();
    let colors: Vec<Color32> = (0..=6).map(|i| figure_color(i, ui)).collect();
    plot.show(ui, |pui| {
        for shape in &fig.shapes {
            match shape {
                Shape::Fill { pts, color } => {
                    let c = colors[(*color).min(6) as usize];
                    let fill = c.gamma_multiply(0.25);
                    if pts.len() <= 5 {
                        pui.polygon(Polygon::new("", PlotPoints::from(pts.clone())).fill_color(fill).stroke(egui::Stroke::NONE));
                    } else {
                        // area under a curve: one strip per sample keeps every piece convex
                        for w in pts[1..pts.len() - 1].windows(2) {
                            let strip = vec![[w[0][0], 0.0], w[0], w[1], [w[1][0], 0.0]];
                            pui.polygon(Polygon::new("", PlotPoints::from(strip)).fill_color(fill).stroke(egui::Stroke::NONE));
                        }
                    }
                }
                Shape::Line { pts, color, dashed, name } => {
                    let mut line = Line::new(name.clone(), PlotPoints::from(pts.clone())).color(colors[(*color).min(6) as usize]).width(if *color >= 6 { 1.0 } else { 2.0 });
                    if *dashed {
                        line = line.style(LineStyle::dashed_loose());
                    }
                    pui.line(line);
                }
                Shape::Arrow { from, to, color, name } => {
                    pui.arrows(Arrows::new(name.clone(), PlotPoints::from(vec![*from]), PlotPoints::from(vec![*to])).color(colors[(*color).min(6) as usize]));
                    pui.line(Line::new(name.clone(), PlotPoints::from(vec![*from, *to])).color(colors[(*color).min(6) as usize]).width(2.0));
                }
                Shape::Dot { at, color, name } => {
                    pui.points(Points::new(name.clone(), PlotPoints::from(vec![*at])).radius(4.5).shape(MarkerShape::Circle).color(colors[(*color).min(6) as usize]));
                }
                Shape::Text { at, text } => {
                    pui.text(Text::new("", PlotPoint::new(at[0], at[1]), RichText::new(text).size(14.0)).color(text_color));
                }
            }
        }
    });
    if !fig.caption.is_empty() {
        ui.label(RichText::new(&fig.caption).small().weak());
    }
}

fn trim(v: f64) -> String {
    let s = format!("{:.6}", v);
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.to_string() }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        egui::Panel::top("modes").show(ui, |ui| {
            ui.add_space(4.0);
            self.mode_bar(ui);
            ui.add_space(4.0);
        });
        match self.tab {
            Tab::Calc => self.calc_tab(ui),
            Tab::Graph => self.graph_tab(ui),
            Tab::Demos => self.demos_tab(ui),
            Tab::Catalog => self.catalog_tab(ui),
        }
    }
}

#[cfg(target_arch = "wasm32")]
fn main() {
    use eframe::wasm_bindgen::JsCast as _;
    wasm_bindgen_futures::spawn_local(async {
        let document = eframe::web_sys::window().and_then(|w| w.document()).expect("no document");
        let canvas = document.get_element_by_id("ticalc_canvas").expect("no canvas").dyn_into::<eframe::web_sys::HtmlCanvasElement>().expect("not a canvas");
        eframe::WebRunner::new().start(canvas, eframe::WebOptions::default(), Box::new(|_cc| Ok(Box::new(App::new())))).await.expect("failed to start");
    });
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    let options = eframe::NativeOptions { viewport: egui::ViewportBuilder::default().with_inner_size([560.0, 780.0]).with_min_inner_size([420.0, 520.0]).with_title("TI-Calc CAS"), ..Default::default() };
    // optional start-up arguments: --eval "expr" (repeatable), --plot "f(x)" (repeatable), --tab graph|catalog
    let mut app = App::new();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut plots = Vec::new();
    for pair in args.chunks(2) {
        match (pair[0].as_str(), pair.get(1)) {
            ("--eval", Some(e)) => {
                app.input = e.clone();
                app.submit();
            }
            ("--plot", Some(f)) => plots.push(PlotFn::new(f)),
            ("--demo", Some(n)) => {
                if let Some(i) = n.parse::<usize>().ok().and_then(|n| n.checked_sub(1)).filter(|i| *i < DEMOS.len()) {
                    app.demo = i;
                    app.demo_steps = cas::demos::run(&DEMOS[i]);
                    app.tab = Tab::Demos;
                }
            }
            ("--tab", Some(t)) => {
                app.tab = match t.as_str() {
                    "graph" => Tab::Graph,
                    "demos" => Tab::Demos,
                    "catalog" => Tab::Catalog,
                    _ => Tab::Calc,
                }
            }
            _ => {}
        }
    }
    if !plots.is_empty() {
        app.fns = plots;
        app.tab = Tab::Graph;
    }
    eframe::run_native("TI-Calc CAS", options, Box::new(|_cc| Ok(Box::new(app))))
}
