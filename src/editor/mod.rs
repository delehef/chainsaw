use anyhow::Context;
use newick::{Newick, NewickTree, NodeID};
use ratatui::{
    backend::Backend,
    crossterm::{
        event::{self, Event, KeyCode, KeyEventKind},
        execute,
        terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
    },
    layout::{Constraint, Direction, Layout},
    prelude::*,
    style::Style,
    widgets::{
        Block, Borders, Cell, Row, Scrollbar, ScrollbarOrientation, ScrollbarState, Table,
        TableState,
    },
    Frame, Terminal,
};
use std::{collections::HashMap, io};

use self::canvas::Canvas;

mod canvas;

enum Screen {
    TreeView,
}

#[derive(Debug)]
struct DispGene {
    name: String,
    species: String,
    landscape: Option<[Vec<String>; 2]>,
}

struct States {
    gene_table: TableState,
    scrollbar: ScrollbarState,
}
impl States {
    fn new(size: usize) -> Self {
        States {
            gene_table: TableState::new().with_selected(0),
            scrollbar: ScrollbarState::new(size - 1),
        }
    }
}

#[derive(Debug)]
enum Clade {
    Taxon {
        graph_line: usize,
        dup_nesting: i16,
        gene: DispGene,
    },
    SubClade {
        subclades: Vec<usize>,
        folded: bool,
    },
}

#[derive(Debug)]
struct CladeHierarchy {
    clades: Vec<Clade>,
}
impl CladeHierarchy {
    fn new() -> Self {
        Self {
            clades: vec![Clade::SubClade {
                subclades: vec![],
                folded: false,
            }],
        }
    }
    fn append_in(&mut self, clade: Clade, parent: usize) -> usize {
        let new = self.clades.len();
        self.clades.push(clade);
        if let Clade::SubClade {
            ref mut subclades, ..
        } = &mut self.clades[parent]
        {
            subclades.push(new);
        } else {
            unreachable!()
        }
        new
    }

    pub fn is_folded(&self, i: usize) -> bool {
        match &self.clades[i] {
            Clade::Taxon { .. } => false,
            Clade::SubClade { folded, .. } => *folded,
        }
    }

    fn get(&self, i: usize) -> &Clade {
        &self.clades[i]
    }

    fn get_mut(&mut self, i: usize) -> &mut Clade {
        &mut self.clades[i]
    }

    fn find_first_taxon(&self, i: usize) -> &Clade {
        match &self.clades[i] {
            taxon @ Clade::Taxon { .. } => taxon,
            Clade::SubClade { subclades, .. } => self.find_first_taxon(subclades[0]),
        }
    }
}

const DEPTH_FACTOR: usize = 2;

struct CanvassedTree {
    canvas: Canvas,
    dup_level: Vec<i16>,
    current_len: usize,
    // screen coordinate -> clade ID
    screen_to_clade: HashMap<usize, Vec<usize>>,
    clades: CladeHierarchy,
}
impl CanvassedTree {
    fn rec_make_tree(
        t: &NewickTree,
        n: NodeID,
        current_y: usize,
        current_depth: usize,
        graph: &mut Canvas,
        dups_nesting: i16,
        dups: &mut Vec<i16>,
        clades: &mut CladeHierarchy,
        current_clade: usize,
    ) -> usize {
        if t[n].is_leaf() {
            dups[current_y] = dups_nesting;
            let my_clade = Clade::Taxon {
                graph_line: current_y,
                dup_nesting: dups_nesting,
                gene: DispGene {
                    name: t.name(n).cloned().unwrap_or("UNKNWN".into()),
                    species: t.attrs(n).get("S").cloned().unwrap_or("UNKNWN".into()),
                    landscape: None,
                },
            };
            clades.append_in(my_clade, current_clade);
            current_y + 1
        } else {
            let my_clade = clades.append_in(
                Clade::SubClade {
                    subclades: vec![],
                    folded: false,
                },
                current_clade,
            );
            let dups_count = if t.is_duplication(n) {
                graph.write_str(current_y, current_depth, "─D");
                dups_nesting + 1
            } else {
                graph.write_str(current_y, current_depth, "─┬");
                dups_nesting
            };

            let old_y = current_y;
            let mut current_y = current_y;
            let bar = if t.is_duplication(n) { '║' } else { '│' };
            let l = if t.is_duplication(n) { '╙' } else { '└' };

            for (i, child) in t[n].children().iter().enumerate() {
                if t[*child].is_leaf() {
                    for i in current_depth + 2..graph.width() {
                        graph.write(current_y, i, '─');
                    }
                }
                for y in old_y + 1..current_y {
                    graph.write(y, current_depth + 1, bar);
                }
                if i > 0 {
                    graph.write(current_y, current_depth + 1, l);
                }

                current_y = Self::rec_make_tree(
                    t,
                    *child,
                    current_y,
                    current_depth + DEPTH_FACTOR,
                    graph,
                    dups_count,
                    dups,
                    clades,
                    my_clade,
                );
            }
            current_y
        }
    }

    fn from_newick(t: &NewickTree) -> Self {
        let leave_count = t.leaves().count();
        let mut canvas = Canvas::new(leave_count, DEPTH_FACTOR * t.topological_depth().1);
        let mut dups = vec![0; leave_count];
        let mut clades: CladeHierarchy = CladeHierarchy::new();

        Self::rec_make_tree(t, t.root(), 0, 0, &mut canvas, 0, &mut dups, &mut clades, 0);

        Self {
            canvas,
            dup_level: dups,
            current_len: leave_count,
            screen_to_clade: Default::default(),
            clades,
        }
    }

    fn len(&self) -> usize {
        self.current_len
    }

    fn to_rows(&mut self) -> Vec<Row> {
        #[derive(Debug)]
        struct RowContext {
            current: usize,
            in_folded: bool,
        }

        let mut rows = Vec::new();
        self.screen_to_clade.clear();

        let mut todos = Vec::new();
        todos.push(RowContext {
            current: 0,
            in_folded: false,
        });

        while let Some(context) = todos.pop() {
            let current = self.clades.get(context.current);
            match current {
                Clade::Taxon {
                    dup_nesting,
                    gene,
                    graph_line,
                } => {
                    self.screen_to_clade.entry(rows.len()).or_default();
                    rows.push(Row::new(vec![
                        Cell::from("│".repeat(*dup_nesting as usize)).light_blue(),
                        self.canvas.line(*graph_line).into(),
                        "".into(),
                        gene.species.clone().into(),
                        gene.name.clone().into(),
                    ]));
                }
                Clade::SubClade { subclades, folded } => {
                    self.screen_to_clade
                        .entry(rows.len())
                        .or_default()
                        .push(context.current);
                    if *folded {
                        if let Clade::Taxon {
                            graph_line,
                            dup_nesting,
                            gene,
                        } = self.clades.find_first_taxon(context.current)
                        {
                            rows.push(Row::new(vec![
                                Cell::from("│".repeat(*dup_nesting as usize)).light_blue(),
                                self.canvas.line(*graph_line).into(),
                                Cell::from("⋮".to_string()).bold().white().on_light_blue(),
                                gene.species.clone().into(),
                                gene.name.clone().into(),
                            ]));
                        }
                    } else {
                        for subclade in subclades.iter().rev() {
                            todos.push(RowContext {
                                current: *subclade,
                                in_folded: context.in_folded,
                            });
                        }
                    }
                }
            }
        }

        self.current_len = rows.len();
        rows
    }

    fn toggle_all_at(&mut self, screen_y: usize) {
        let target_state = !self.screen_to_clade[&screen_y]
            .iter()
            .any(|c| self.clades.is_folded(*c));

        for clade in self.screen_to_clade.get(&screen_y).unwrap().iter().rev() {
            if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
                *folded = target_state;
            } else {
                unreachable!()
            }
        }
    }

    fn fold_at(&mut self, screen_y: usize) {
        for clade in self.screen_to_clade.get(&screen_y).unwrap().iter().rev() {
            if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
                if !*folded {
                    *folded = true;
                    return;
                }
            } else {
                unreachable!()
            }
        }
    }

    fn unfold_at(&mut self, screen_y: usize) {
        for clade in self.screen_to_clade.get(&screen_y).unwrap().iter() {
            if let Clade::SubClade { ref mut folded, .. } = self.clades.get_mut(*clade) {
                if *folded {
                    *folded = false;
                    return;
                }
            } else {
                unreachable!()
            }
        }
    }

    fn max_dup_nesting(&self) -> u16 {
        self.dup_level.iter().max().cloned().unwrap_or(0) as u16
    }
}

struct Editor {
    name: String,
    tree: NewickTree,
    plot: CanvassedTree,
    screen: Screen,
    states: States,
}
impl Editor {
    pub fn new(name: String, tree: NewickTree) -> Self {
        let plot = CanvassedTree::from_newick(&tree);
        let states = States::new(plot.len());

        Self {
            name,
            tree,
            plot,
            screen: Screen::TreeView,
            states,
        }
    }

    fn move_to(&mut self, i: usize) {
        self.states.gene_table.select(Some(i));
        self.states.scrollbar = self.states.scrollbar.position(i);
    }

    fn prev(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| i.saturating_sub(count))
            .unwrap_or_default();
        self.move_to(i);
    }

    fn next(&mut self, count: usize) {
        let i = self
            .states
            .gene_table
            .selected()
            .map(|i| (i + count).clamp(0, self.plot.len() - 1))
            .unwrap_or_default();
        self.move_to(i);
    }

    fn top(&mut self) {
        self.move_to(0);
    }

    fn bottom(&mut self) {
        self.move_to(self.plot.len() - 1);
    }

    fn render(&mut self, f: &mut Frame) {
        const INDENT_FACTOR: usize = 2;
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Min(1),
                Constraint::Length(3),
            ])
            .split(f.size());

        let title_block = Block::default()
            .borders(Borders::BOTTOM)
            .style(Style::default())
            .title(Line::from(vec![
                "←".yellow().bold(),
                " fold 1×  ".into(),
                "→".yellow().bold(),
                " unfold 1×  ".into(),
                "[TAB]".yellow().bold(),
                " cycle fold  ".into(),
            ]));

        // let title = Paragraph::new(Text::styled(&self.name, Style::default())).block(title_block);
        f.render_widget(title_block, chunks[0]);

        let tree_depth = self.tree.topological_depth().1;
        let dups_width = self.plot.max_dup_nesting();
        let rows = self.plot.to_rows();
        let widths = [
            Constraint::Length(dups_width),
            Constraint::Length((INDENT_FACTOR * tree_depth) as u16),
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Fill(1),
        ];
        let table = Table::new(rows, widths)
            .column_spacing(1)
            .header(
                Row::new(vec!["Dup.", "", "", "Species", "Gene"])
                    .style(Style::new().bold())
                    .bottom_margin(1),
            )
            // The selected row and its content can also be styled.
            .highlight_style(Style::new().reversed());

        f.render_stateful_widget(table, chunks[1], &mut self.states.gene_table);

        f.render_stateful_widget(
            Scrollbar::default()
                .orientation(ScrollbarOrientation::VerticalRight)
                .begin_symbol(Some("^"))
                .end_symbol(Some("v")),
            chunks[1].inner(Margin {
                vertical: 1,
                horizontal: 1,
            }),
            &mut self.states.scrollbar,
        );
    }

    fn run<B: Backend>(&mut self, terminal: &mut Terminal<B>) -> anyhow::Result<()> {
        loop {
            terminal.draw(|term| self.render(term))?;
            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    match key.code {
                        KeyCode::Char('q') => return Ok(()),
                        KeyCode::Up => self.prev(1),
                        KeyCode::Down => self.next(1),
                        KeyCode::PageUp => self.prev(10),
                        KeyCode::PageDown => self.next(10),
                        KeyCode::Home => self.top(),
                        KeyCode::End => self.bottom(),
                        KeyCode::Left => {
                            self.plot
                                .fold_at(self.states.gene_table.selected().unwrap());
                        }
                        KeyCode::Right => {
                            self.plot
                                .unfold_at(self.states.gene_table.selected().unwrap());
                        }
                        KeyCode::Tab => self
                            .plot
                            .toggle_all_at(self.states.gene_table.selected().unwrap()),
                        _ => {}
                    }
                    terminal.draw(|term| self.render(term))?;
                }
            }
        }
    }
}

pub fn run(name: String, t: NewickTree) -> anyhow::Result<()> {
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen).context("failed to initialize terminal")?;
    let backend = ratatui::backend::CrosstermBackend::new(stdout);
    let mut terminal = ratatui::Terminal::new(backend)?;

    let mut editor = Editor::new(name, t);
    editor.run(&mut terminal)?;

    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen,)?;
    terminal.show_cursor()?;

    Ok(())
}
