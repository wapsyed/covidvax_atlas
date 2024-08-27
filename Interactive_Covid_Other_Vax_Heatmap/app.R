
load('htShiny.RData')
chooseCRANmirror(ind = 1)
setRepositories(ind = 1:4)
if(!requireNamespace('InteractiveComplexHeatmap', quietly = TRUE)) {
	install.packages('InteractiveComplexHeatmap')
}
library(shiny)
suppressPackageStartupMessages(library(InteractiveComplexHeatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))
ui = fluidPage(
  h1("COVID-19 Vaccination Atlas - Gene explorer"),
  h2("Comparison with other vaccines - Immune genes"),
  p("Reference: COVID-19 Vaccination Atlas: An Integrative Systems Vaccinology Approach"),
  a(href = "https://github.com/wapsyed/covidvax_atlas", "Access the GitHub Repository"),
  p("Explore immune genes by either selecting an area or typing a list of genes of interest (e.g. 'IL2, IL3, IL10...') in the 'magnifier' icon below the left panel."),
  p("To find a family of genes that show a pattern in their names, such as interleukins (IL*), select the 'Regular expression' box and type '^IL' to show genes that start with IL."),
  p("Otherwise, it will show all the genes that have 'IL' in their names."),
  p("The left annotations correspond to the gene role in innate (shades of green) and adaptive immune (shades of purple) systems, as well as other general immune processes, such as the complement system."),
  p("To visualize all gene annotations, select the genes by either of the options above, click on the table icon under the 'Selected sub-heatmap', and select 'open table'."),
  
  hr(), 
  InteractiveComplexHeatmapOutput(
    heatmap_id = heatmap_id, 
    title1 = title1, 
    title2 = title2,
    width1 = 700, 
    height1 = 1000, 
    width2 = 600, 
    height2 = height2, 
    layout = layout, 
    compact = compact,
    action = action, 
    cursor = cursor, 
    response = response, 
    brush_opt = brush_opt, 
    output_ui_float = output_ui_float
  ), 
	html
)

server = function(input, output, session) {
	makeInteractiveComplexHeatmap(input, output, session, ht_list, 
		show_cell_fun = show_cell_fun, show_layer_fun = show_layer_fun)
}

shinyApp(ui = ui, server = server)
