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
                  tabPanel("Overview",
                           HTML("
 <h3>Welcome to ProGeny SSP Generator!</h3>
 <p>
 The <strong>ProGeny simple stellar population (SSP) generator</strong> can be used to create highly flexible SSPs by incorporating the user's choice of isochrones, stellar spectral atmospheres, and IMFs. The various tabs reflect the underlying stages of SSP generation and give users access to nearly 100% of ProGeny's functionality—without needing to learn specific R syntax.
 </p>
 <p>
 The goal is to allow users to generate SSPs and use them directly in ProSpect (as FITS files), or process them for your SED code of choice. For deeper understanding, please refer to the associated papers: Robotham & Bellstedt (2025, RASTI) and Bellstedt & Robotham (2025, MNRAS).
 </p>

 <p>
 The sub sections below provide basic information on the related tabs, and in general users should progress their ProGeny journey from left to right (so Isochrone -> Atmospheres -> Interpolate -> Make SSP). The additional tabs beyond that provide a bit more information about the data underlying the Isochrone and Atmospheres tabs. There is much more information available through the ProGeny papers and R package itself, so if you need much more information that this you should look there (before getting in touch with questions those cover).
</p>

<p>
 If you have specific questions about ProGeny that is not well covered in the papers and the R package documentation please get in touch at <a href='mailto:aaron.robotham@uwa.edu.au'>aaron.robotham@uwa.edu.au</a>.
 </p>

 <h4>Isochrone</h4>
 <p>
 Isochrones encode how stars evolve over cosmic time. Users should select one from the dropdown list. The <em>Iso Info</em> tab provides more details. We recommend starting here. Our preference is the MIST isochrones due to their broad coverage, but experimentation is encouraged.
 </p>

 <p>
 The selected isochrone should load pretty quickly (a few seconds at worst), and you should see a useful Temp versus logG figure to give you an idea of the evolutionary coverage of the chosen isochrone.
 </p>

 <h4>Atmospheres</h4>
 <p>
 Stellar spectral atmospheres provide the spectral data for stars across different physical properties. The <em>Atmos Info</em> tab offers more detail. Atmosphere types include:
 <ul>
 <li><strong>base</strong> and <strong>extend</strong>: main sequence and some evolved phases</li>
 <li><strong>AGB</strong>: asymptotic giant branch stars</li>
 <li><strong>hot</strong>: young and horizontal branch stars</li>
 <li><strong>white</strong>: white dwarfs</li>
 <li><strong>WR</strong>: Wolf-Rayet stars</li>
 </ul>
 Users can upload their own spectral grid (a single-column text file of wavelengths in Angstroms, no header). The default is the high-resolution CB19 gridding, which balances resolution and performance.
 </p>

 <p>
 The selected atmospheres should take a few tens of seconds to load, and you should see a couple of userful plots. One will overlay and colour-code the various atmosphere grids on the isochrone Temp versus logG figure to give you an idea of the evolutionary coverage of the chosen isochrone and how well sampled the properties are with your chosen stellar atmospheres. The other shows the spectral resolution as a functions of wavelength. If users do not upload their own this will simply be the default high-resolution CB19 grid.
 </p>

 <h4>Interpolate</h4>
 <p>
 This tab controls how stellar properties are interpolated from the gridded atmospheres onto the non-gridded isochrones. The defaults are robust and were used in the reference papers, but users can experiment to explore their effects. Almost certainly the auto options regarding the labels should be left as they are, unless you really know what you are doing (and how these work).
 </p>

 <p>
 The interpolation should take a few tens of seconds to run There is nothing pretty generated by this tab, rather it will just tell you it has succesfully run.
 </p>

 <h4>IMF</h4>
 <p>
 The initial mass function (IMF) determines the distribution of stellar masses. ProGeny uniquely supports evolving IMFs. The default is a Chabrier IMF, but this may not suit all science cases (e.g., high redshift). This is a critical tab—users should carefully consider their choice here. Note how the available options vary depending on the IMF selected? That is partly me pointing it out because I think it is very neat, but it is also relevant because some of the parameters have different effects depending on the IMF. Check out the auto-generated plot to better understand the impact of your options on the resulting IMF.
 </p>

 <h4>Make SSP</h4>
 <p>
 This tab runs the SSP generation, combining outputs from the other tabs. Users should make sure all the previous tabs have been run, and the assoicated text on this tab is TRUE (helpfully also coloured green). It may take a minute or two to complete the main SSP building even on 8 cores (so please be patient).
 </p>

 <p>
 After building your SSP, users can run 'Check SSP' to run lots of internally checks to ensure the SSP seems to be valid and legal for running in ProSpect. This is the minimum grade for checking the SSP is appropriate for science, but really it is only checking the data format looks sensible. If users have made poor choices when building the various input for ProGeny (the previous tabs) than there might be scientific issues still.
 </p>

 <p>
 Once you have run your checks, the SSP is ready to be downloaded with the 'Download SSP' button. The output will be FITS files in the standard ProSpect SSP format. This is detailed fully in the ProSpect package, but in brief see below for the included BC03lr data. Within R you can load a generated ProGeny SSP with the ProSpect convenience function speclib_FITSload, but in general the file is a pretty simple multi-extension FITS that can be loaded directly into any software that supports that format (including, of course, AstroPy etc).
 </P>

<h5>Overview: Low Resolution BC03lr Format</h5>
 The highest-level list structure contains the following elements:
 <ul>
 <li><strong>Z</strong>: Numeric vector of length 6 describing the availability and location of different metallicity stellar populations.</li>
 <li><strong>Age</strong>: Numeric vector of length 221 giving the midpoint of the age in years for the stellar populations.</li>
 <li><strong>AgeBins</strong>: Numeric vector of length 222 giving the age bin extremes in years of the stellar populations.</li>
 <li><strong>AgeWeights</strong>: Numeric vector of length 221 giving the total time in years for the stellar population bins.</li>
 <li><strong>Wave</strong>: Numeric vector of length 1221 giving the spectral wavelength in Angstroms.</li>
 <li><strong>Labels</strong>: Character vector of length 5 giving handy axis labels.</li>
 <li><strong>Zspec</strong>: List of 6 (for each metallicity) 221x1221 matrices that give the full spectral information.</li>
 <li><strong>Zevo</strong>: List of 6 (for each metallicity) 221x5 data frames that contain the evolution information for the stellar material:
 <ul>
 <li>SMstar: Stellar mass in stars</li>
 <li>SMgas: Stellar mass in gas</li>
 <li>SMtot: Total stellar mass</li>
 <li>SFR: Star formation rate</li>
 <li>SMrem: Stellar mass in remnants</li>
 </ul>
 </li>
 </ul>
 ")
                  ),
                  tabPanel("Isochrone",
                           verbatimTextOutput("iso_status"),
                           plotOutput("plot_iso", height = "600px"),
                           verbatimTextOutput("iso_summary")
                  ),
                  tabPanel("Atmospheres",
                           verbatimTextOutput("selectedPath"),
                           verbatimTextOutput("atmos_status"),
                           plotOutput("plot_atmos", height = "600px"),
                           #verbatimTextOutput("loaded_wave_samp"),
                           tags$h4("Wavelength grid info below:"),
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
