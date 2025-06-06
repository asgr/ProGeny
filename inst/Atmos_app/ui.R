ui = fluidPage(
  shinybusy::use_busy_spinner(spin = "fading-circle"),
  titlePanel("ProGeny SSP Generator"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tabs == 'Isochrone'",
        tags$h5("Allow 10 sec to load Iso File"),
        fileInput("iso_file", "Choose Isochrone File [.fst]", accept = ".fst"),
        br(),
        actionButton("iso_done", "Return Isochrone"),
      ),

      conditionalPanel(
        condition = "input.tabs == 'Atmospheres'",
        shinyFiles::shinyDirButton("destpath", "Atmos Path", "Select a folder"),
        selectInput("base", "Base", choices = c("C3K Conroy" = "combine_C3K_conroy", "PHOENIX Husser" = "combine_PHOENIX_husser", "PHOENIX Allard" = "combine_PHOENIX_Allard", "MILES Vazdekis" = "combine_MILES_vazdekis", "BaSeL-3.1 WLBC" = "combine_BASEL_wlbc")),
        selectInput("extend", "Extend", choices = c("PHOENIX Allard" = "combine_PHOENIX_Allard", "ATLAS9 Castelli" = "combine_ATLAS9_castelli", "None" = "None")),
        selectInput("hot", "Hot", choices = c("PoWR" = "combine_OB_PoWR", "hot" = "combine_hot", "None" = "None")),
        selectInput("AGB", "AGB", choices = c("AGB Lancon" = "combine_AGB_lancon", "None" = "None")),
        selectInput("white", "White Dwarfs", choices = c("TMAP Werner" = "combine_TMAP_werner", "white" = "combine_white", "None" = "None")),
        selectInput("WR", "Wolf-Rayet", choices = c("PoWR" = "combine_WNE_PoWR", "None" = "None")),
        fileInput("wave_file", "User Wave (Ang) [.tab .dat .txt]", accept = c('tab', 'dat', 'txt')),
        numericInput("atmos_cores", "Number of Cores", value = 8, min = 1, step=1),
        tags$h5("Allow 30 sec to [Load Atmos]"),
        actionButton("load_atmos", "Load Atmos"),
        br(), br(),
        actionButton("atmos_done", "Return Atmos")
      ),

      conditionalPanel(
        condition = "input.tabs == 'Interpolate'",
        sliderInput('radius', 'Search Radius', min=1, max=4, value=2),
        sliderInput('weight_pow', 'Weight Power', min=1, max=4, value=2),
        sliderInput('k', 'Search k', min=4, max=16, value=8, step=1),
        checkboxInput('do_hot', 'Use Hot', value = TRUE),
        checkboxInput('do_AGB', 'Use AGB', value = TRUE),
        checkboxInput('do_white', 'Use White Dwarf', value = TRUE),
        checkboxInput('do_WR', 'Use Wolf-Rayet', value = TRUE),
        selectInput('label_AGB', 'AGB Phase', choices=-1:10, multiple=TRUE),
        selectInput('label_white', 'White Dwarf Phase', choices=-1:10, multiple=TRUE),
        selectInput('label_WR', 'AGB Wolf-Rayet', choices=-1:10, multiple=TRUE),
        tags$h5("Allow 30 sec to [Run Interp]"),
        actionButton("run_interp", "Run Interp"),
        br(), br(),
        actionButton("interp_done", "Return Interp Grids"),
      ),

      conditionalPanel(
        condition = "input.tabs == 'IMF'",
        selectInput("imf", "IMF", choices = c("Chabrier" = "IMF_Chabrier",
                                              "Kroupa" = "IMF_Kroupa",
                                              "Salpeter" = "IMF_Salpeter",
                                              "Kroupa_evo" = "IMF_Kroupa_evo",
                                              "Kroupa_Zevo" = "IMF_Kroupa_Zevo",
                                              "Lacey_evo" = "IMF_Lacey_evo",
                                              "Lacey_Zevo" = "IMF_Lacey_Zevo")),
        tags$h4("IMF parameters:"),
        uiOutput("dynamic_imf")
      ),

      conditionalPanel(
        condition = "input.tabs == 'Make SSP'",
        numericInput("SSP_cores", "Number of Cores", value = 8, min = 1, step=1),
        tags$h5("Allow 2 min to [Make SSP]"),
        actionButton("make_ssp", "Make SSP"),
        br(), br(),
        actionButton("check_ssp", "Check SSP"),
        br(), br(),
        actionButton("return_ssp", "Return SSP"),
        br(), br(),
        downloadButton("download_ssp", "Download SSP [.fits]")
      )
    ),

    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Isochrone",
                           verbatimTextOutput("iso_status"),
                           plotOutput("plot_iso", height = "600px"),
                           verbatimTextOutput("iso_summary")
                  ),
                  tabPanel("Atmospheres",
                           verbatimTextOutput("selectedPath"),
                           verbatimTextOutput("loaded_wave_samp"),
                           verbatimTextOutput("atmos_status"),
                           plotOutput("plot_atmos", height = "600px"),
                           verbatimTextOutput("summary_wave_samp"),
                           plotOutput("plot_wave_samp", height = "400px"),
                  ),
                  tabPanel("Interpolate",
                           verbatimTextOutput("interp_status")
                  ),
                  tabPanel("IMF",
                           plotOutput("plot_imf", height = "600px")
                  ),
                  tabPanel("Make SSP",
                           tags$h4("Pre-requisite Inputs Status:"),
                           uiOutput("used_imf"),
                           tags$h4("[Make SSP] Status:"),
                           uiOutput("SSP_status"),
                           tags$h4("[Check SSP] Status:"),
                           uiOutput("SSP_check"),
                           tags$h4("[Return SSP] Status:"),
                           uiOutput("SSP_return"),
                           tags$h4("[Download SSP] Status:"),
                           uiOutput("SSP_return"),
                  ),

                  tabPanel("Iso Info",
                           HTML("
<h3>Available ProGeny Isochrone Libraries</h3>
 <p>All libraries are in <strong>FST</strong> format:</p>
 <ul>
 <li><strong>ParsecIso</strong> [fst]</li>
 <li><strong>MistIso</strong> [fst]</li>
 <li><strong>BastiIso_New</strong> [fst]</li>
 <li><strong>BastiIso_FSPS</strong> [fst]</li>
 </ul>

 <h4><strong>MistIso</strong></h4>
 <p>The Mist isochrones were obtained from
 <a href='https://waps.cfa.harvard.edu/MIST/model_grids.html' target='_blank'>
 MIST model grids</a>. We use the pre-packaged isochrone with v/v<sub>crit</sub> = 0.4,
 and downloaded the basic version of the catalogue (25 columns).</p>
 <ul>
 <li><strong>MIST version:</strong> 1.2</li>
 <li><strong>MESA revision:</strong> 7503</li>
 </ul>
 <pre>
---------------------------------------------
 Yinit        Zinit   [Fe/H]  [a/Fe]  v/vcrit
0.2490  1.42857E-06.   -4.00.   0.00.    0.40
---------------------------------------------
 </pre>
 <p><strong>Number of isochrones:</strong> 107</p>
 <p><em>General comment:</em> Mist tracks remnants, but they appear to be missing for older, low-metallicity isochrones. Likely usable range is logZ ≥ -2. Covers 15 Z values from 0.000002 to 0.06 (x3 solar) — the broadest Z range.</p>

 <hr>

 <h4><strong>ParsecIso</strong></h4>
 <p>Generated by CMD 3.7 (<a href='http://stev.oapd.inaf.it/cmd' target='_blank'>link</a>) on 2024-06-07. Based on PARSEC v1.2S + COLIBRI S_37/S_35/PR16. Includes thermal pulse cycles, Reimers mass loss (η=0.2), LPV periods, and YBC bolometric corrections.</p>
 <p>Photometric system: <i>UBVRIJHK</i>. Includes circumstellar dust models (O-rich and C-rich).</p>
 <p><em>General comment:</em> Missing remnants almost entirely. No UV upturn for old populations. 27 Z values from 0.0001 to 0.06 — the most densely Z-sampled library.</p>

 <hr>

 <h4><strong>BastiIso_New</strong></h4>
 <p>From the <a href='http://basti-iac.oa-abruzzo.inaf.it/isocs.html' target='_blank'>BaSTI-IAC database</a>. Downloaded 2024-09-12. logG computed from Lum, Mass, and Teff using Stefan-Boltzmann Law (see <code>BaSTI_New_isochrone.R</code>).</p>
 <p><em>General comment:</em> Missing remnants. No UV upturn for old populations. 10 Z values from 0.0003 to 0.04 — narrowest and least sampled Z range.</p>

 <hr>

 <h4><strong>BastiIso_FSPS</strong></h4>
 <p>Converted BaSTI isochrones used in FSPS. logG values added manually for remnants (originally set to -99). Based on older BaSTI versions (~2010–2015).</p>
 <p><em>General comment:</em> Remnants added by hand. Provides reasonable UV upturn. 21 Z values from 0.00001 to 0.06 — middling in Z coverage and sampling.</p>
")
                  ),

                  tabPanel("Atmos Info",
                           HTML("
 <h3>Available ProGeny Atmosphere Libraries</h3>
 <p>All libraries are in <strong>FITS</strong> format:</p>

 <h4>Base</h4>
 <ul>
 <li><strong>combine_C3K_conroy</strong> [preferred]</li>
 <li>combine_PHOENIX_husser</li>
 <li>combine_PHOENIX_husser_4MOST</li>
 <li>combine_PHOENIX_husser_NIRspec</li>
 <li>combine_MILES_vazdekis</li>
 <li>combine_BASEL_wlbc</li>
 </ul>

 <h4>Extend</h4>
 <ul>
 <li><strong>combine_PHOENIX_allard</strong> [preferred]</li>
 <li>combine_PHOENIX_allard_4MOST</li>
 <li>combine_PHOENIX_allard_NIRspec</li>
 <li>combine_ATLAS9_castelli</li>
 </ul>

 <h4>Hot</h4>
 <ul>
 <li><strong>combine_OB_PoWR</strong> [preferred]</li>
 <li>combine_hot <em>(not for science)</em></li>
 </ul>

 <h4>AGB</h4>
 <ul>
 <li><strong>combine_AGB_lancon</strong> [preferred]</li>
 </ul>

 <h4>White Dwarfs</h4>
 <ul>
 <li><strong>combine_TMAP_werner</strong> [preferred]</li>
 <li>combine_white</li>
 </ul>

 <h4>Wolf-Rayet</h4>
 <ul>
 <li><strong>combine_WNE_PoWR</strong> [preferred 'WR']</li>
 </ul>

 <hr>

 <h4>General Notes</h4>
 <ul>
 <li>All libraries normalized so ∫f(λ)dλ = 1 — scale by luminosity for isochrone spectra.</li>
 <li>Re-binned to CB19 resolution (unless marked _4MOST or _NIRspec): 90 Å – 3.6 cm, median 0.5 Å, 16902 λ-points.</li>
 <li>Native wavelength ranges and sampling vary — see individual library documentation.</li>
 <li>Other libraries (e.g. Levenhagen 2017) excluded due to insufficient wavelength coverage.</li>
 </ul>
 ")

                  )
      )
    )
  )
)
