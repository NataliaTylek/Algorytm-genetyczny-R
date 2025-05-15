library(osrm)
library(dplyr)
library(maps)
library(ggplot2)
library(grid)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


poland_sf <- ne_countries(scale = "medium", country = "Poland", returnclass = "sf")

data(world.cities)

polskie_miasta <- world.cities %>%
  filter(country.etc == "Poland") %>%
  arrange(desc(pop))

create_distance_matrix <- function(miasta) {
  miasta_osrm <- miasta %>% select(long, lat) %>% as.data.frame()
  rownames(miasta_osrm) <- miasta$id
  macierz <- osrmTable(loc = miasta_osrm)
  dystanse <- macierz$durations
  return(dystanse)
}

miasta_16 <- head(polskie_miasta, 16) %>% select(name, lat, long)
miasta_30 <- head(polskie_miasta, 30) %>% select(name, lat, long)
miasta_50 <- head(polskie_miasta, 50) %>% select(name, lat, long)
miasta_80 <- head(polskie_miasta, 80) %>% select(name, lat, long)

miasta_16$id <- miasta_16$name
miasta_30$id <- miasta_30$name
miasta_50$id <- miasta_50$name
miasta_80$id <- miasta_80$name

macierz_16 <- create_distance_matrix(miasta_16)
macierz_30 <- create_distance_matrix(miasta_30)
macierz_50 <- create_distance_matrix(miasta_50)
macierz_80 <- create_distance_matrix(miasta_80)


plot_city_map <- function(miasta, title) {
  ggplot() +
    geom_sf(data = poland_sf, fill = NA, color = "black") +
    geom_point(data = miasta, aes(x = long, y = lat), color = "blue", size = 3) +
    geom_text(data = miasta, aes(x = long, y = lat, label = id), vjust = -0.5, size = 3) +
    labs(title = title, x = "Długość geograficzna", y = "Szerokość geograficzna") +
    coord_sf() +
    theme_minimal() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))
}

plot_city_map(miasta_16, "Rozmieszczenie 16 największych miast Polski")
plot_city_map(miasta_30, "Rozmieszczenie 30 największych miast Polski")
plot_city_map(miasta_50, "Rozmieszczenie 50 największych miast Polski")
plot_city_map(miasta_80, "Rozmieszczenie 80 największych miast Polski")


plot_best_route_on_map <- function(best_route, miasta, title) {
  route_coords <- miasta[best_route, c("long", "lat", "name")]
  route_coords$order <- 1:nrow(route_coords)
  route_coords <- rbind(route_coords, route_coords[1, ])
  route_coords$order[nrow(route_coords)] <- nrow(route_coords)
  
  ggplot() +
    geom_sf(data = poland_sf, fill = NA, color = "black") +
    geom_point(data = miasta, aes(x = long, y = lat), color = "blue", size = 3) +
    geom_point(data = route_coords[1, ], aes(x = long, y = lat), color = "green", size = 4) +
    geom_text(data = route_coords[-nrow(route_coords), ],
              aes(x = long, y = lat, label = order), vjust = -0.7, size = 3.5, color = "black") +
    geom_text(data = miasta, aes(x = long, y = lat, label = name), vjust = 1.5, size = 2.8, color = "gray40") +
    geom_path(data = route_coords, aes(x = long, y = lat),
              color = "red", linetype = "dashed", size = 1,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"))) +
    labs(title = title, x = "Długość geograficzna", y = "Szerokość geograficzna") +
    theme_minimal() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))
}

calculate_distance <- function(route, distance_matrix) {
  total_distance <- sum(sapply(1:(length(route) - 1), function(i) distance_matrix[route[i], route[i+1]]))
  total_distance <- total_distance + distance_matrix[route[length(route)], route[1]]
  return(total_distance)
}

genetic_algorithm <- function(distance_matrix,
                              num_generations = 400,
                              pop_size = 80,
                              mutation_prob = 0.1,
                              elite_size = 2,
                              tournament_size = 3) {
  
  num_cities <- nrow(distance_matrix)
  fitness_history <- numeric(num_generations)
  current_history <- numeric(num_generations)
  
  crossover <- function(parent1, parent2) {
    pt <- sort(sample(1:length(parent1), 2))
    child <- rep(NA, length(parent1))
    child[pt[1]:pt[2]] <- parent1[pt[1]:pt[2]]
    remaining <- parent2[!parent2 %in% child]
    child[is.na(child)] <- remaining
    return(child)
  }
  
  mutate <- function(route) {
    idx <- sample(1:length(route), 2)
    route[idx] <- route[rev(idx)]
    return(route)
  }
  
  tournament_selection <- function(population, distances, tournament_size) {
    candidates_idx <- sample(1:length(population), tournament_size)
    best_idx <- candidates_idx[which.min(distances[candidates_idx])]
    return(population[[best_idx]])
  }
  
  population <- replicate(pop_size, sample(1:num_cities), simplify = FALSE)
  best_solution <- NULL
  best_distance <- Inf
  
  for (generation in 1:num_generations) {
    distances <- sapply(population, calculate_distance, distance_matrix = distance_matrix)

    elite_indices <- order(distances)[1:elite_size]
    elite_individuals <- population[elite_indices]
    
    new_population <- list()
    
    for (i in seq(elite_size + 1, pop_size, by = 2)) {
      p1 <- tournament_selection(population, distances, tournament_size)
      p2 <- tournament_selection(population, distances, tournament_size)
      c1 <- crossover(p1, p2)
      c2 <- crossover(p2, p1)
      if (runif(1) < mutation_prob) c1 <- mutate(c1)
      if (runif(1) < mutation_prob) c2 <- mutate(c2)
      new_population[[i]] <- c1
      if (i + 1 <= pop_size) new_population[[i + 1]] <- c2
    }
    
    for (i in 1:elite_size) {
      new_population[[i]] <- elite_individuals[[i]]
    }

    population <- new_population
    
    distances <- sapply(population, calculate_distance, distance_matrix = distance_matrix)
    gen_best_distance <- min(distances)
    fitness_history[generation] <- gen_best_distance
    current_history[generation] <- mean(distances)
    
    if (gen_best_distance < best_distance) {
      best_distance <- gen_best_distance
      best_solution <- population[[which.min(distances)]]
    }
  }
  
  return(list(
    best_solution = best_solution,
    best_distance = best_distance,
    fitness_history = fitness_history,
    current_history = current_history
  ))
}



start_time <- Sys.time()

result_ga <- genetic_algorithm(macierz)

end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))


ggplot(data.frame(Iteracja = 1:length(result_ga$fitness_history),
                  Długość = result_ga$fitness_history),
       aes(x = Iteracja, y = Długość)) +
  geom_line(color = "blue") +
  labs(title = "Postęp GA – długość trasy vs iteracja", x = "Iteracja", y = "Długość trasy") +
  theme_minimal()


  rozmiary <- c(16, 30, 50, 80)
kolory <- list(
  "16" = c("darkblue", "lightblue"),
  "30" = c("darkred", "lightcoral"),
  "50" = c("darkgreen", "lightgreen"),
  "80" = c("darkorange", "moccasin")
)

wyniki_ga <- list()
df_list_ga <- list()
mapy_ga <- list()

for (rozmiar in rozmiary) {
  macierz <- get(paste0("macierz_", rozmiar))
  miasta <- get(paste0("miasta_", rozmiar))
  
  start_time <- Sys.time()
  result_ga <- genetic_algorithm(macierz)
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  wyniki_ga[[as.character(rozmiar)]] <- result_ga
  
  cat(sprintf("Rozmiar: %d\nNajlepsza trasa (GA): %s\nNajlepsza długość trasy: %.2f\nCzas działania: %.3f sekund\n\n",
              rozmiar,
              paste(result_ga$best_solution, collapse = "-"),
              result_ga$best_distance,
              execution_time))
  
  mapy_ga[[as.character(rozmiar)]] <- plot_best_route_on_map(
    result_ga$best_solution,
    miasta,
    paste0("Najlepsza trasa - Algorytm Genetyczny (", rozmiar, " miast)")
  )
  
  n <- which(result_ga$fitness_history == 0)[1]
  if (is.na(n)) n <- length(result_ga$fitness_history)
  
  df_temp <- data.frame(
    Iteracja = 1:(n),
    Najlepsza = result_ga$fitness_history[1:n],
    Aktualna = result_ga$current_history[1:n],
    Rozmiar = as.factor(rozmiar)
  )
  
  df_list_ga[[as.character(rozmiar)]] <- df_temp
}
wykresy_ga <- lapply(rozmiary, function(rozmiar) {
  df_temp <- df_list_ga[[as.character(rozmiar)]]
  kolory_rozmiaru <- kolory[[as.character(rozmiar)]]
  
  ggplot(df_temp, aes(x = Iteracja)) +
    geom_line(aes(y = Aktualna), color = kolory_rozmiaru[2]) +
    geom_line(aes(y = Najlepsza), color = kolory_rozmiaru[1]) +
    labs(title = paste0("Postęp GA – ", rozmiar, " miast"),
         x = "Iteracja", y = "Długość trasy") +
    theme_minimal()
})
do.call(grid.arrange, c(wykresy_ga, ncol = 1))

print(mapy_ga[["16"]])
print(mapy_ga[["30"]])
print(mapy_ga[["50"]])
print(mapy_ga[["80"]])


