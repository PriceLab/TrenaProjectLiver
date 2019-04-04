library(TrenaViz)
tv <- TrenaViz("TrenaProjectLiver")
runApp(createApp(tv, port=6838))
later(function(){browseURL("http://0.0.0.0:6838")}, 2)
