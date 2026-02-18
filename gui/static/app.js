/**
 * RealDetect Web GUI — Map Application
 *
 * Leaflet-based tiled map with real-time station and event display.
 * Communicates with the FastAPI backend over WebSocket.
 */

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let map;
let stationLayerGroup;
let eventLayerGroup;
let arrivalLayerGroup;
let stations = [];
let events = [];
let selectedEventId = null;
let ws;
let wsReconnectTimer = null;

// Base map tile layers
const baseLayers = {
  osm: {
    url: "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
    attr: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>',
  },
  topo: {
    url: "https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png",
    attr: '&copy; <a href="https://opentopomap.org">OpenTopoMap</a>',
  },
  dark: {
    url: "https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png",
    attr: '&copy; <a href="https://carto.com/">CARTO</a>',
  },
  satellite: {
    url: "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    attr: '&copy; <a href="https://www.esri.com/">Esri</a>',
  },
};

let currentBaseLayer = null;

// ---------------------------------------------------------------------------
// Initialize map
// ---------------------------------------------------------------------------
function initMap() {
  map = L.map("map", {
    center: [35, -100],
    zoom: 4,
    zoomControl: true,
    preferCanvas: true,
  });

  setBaseLayer("dark");

  stationLayerGroup = L.layerGroup().addTo(map);
  eventLayerGroup = L.layerGroup().addTo(map);
  arrivalLayerGroup = L.layerGroup();

  // Layer toggle listeners
  document.getElementById("layer-stations").addEventListener("change", function () {
    this.checked ? stationLayerGroup.addTo(map) : stationLayerGroup.remove();
  });
  document.getElementById("layer-events").addEventListener("change", function () {
    this.checked ? eventLayerGroup.addTo(map) : eventLayerGroup.remove();
  });
  document.getElementById("layer-arrivals").addEventListener("change", function () {
    this.checked ? arrivalLayerGroup.addTo(map) : arrivalLayerGroup.remove();
  });
  document.getElementById("basemap-select").addEventListener("change", function () {
    setBaseLayer(this.value);
  });
}

function setBaseLayer(name) {
  if (currentBaseLayer) map.removeLayer(currentBaseLayer);
  const cfg = baseLayers[name] || baseLayers.dark;
  currentBaseLayer = L.tileLayer(cfg.url, {
    attribution: cfg.attr,
    maxZoom: 18,
  }).addTo(map);
}

// ---------------------------------------------------------------------------
// Station rendering
// ---------------------------------------------------------------------------
function renderStations(data) {
  stations = data;
  stationLayerGroup.clearLayers();

  stations.forEach(function (st) {
    const marker = L.circleMarker([st.latitude, st.longitude], {
      radius: 5,
      fillColor: "#4e9af1",
      fillOpacity: 0.85,
      color: "#7ab8ff",
      weight: 1,
    });

    marker.bindTooltip(
      '<div class="station-tooltip"><b>' +
        st.network +
        "." +
        st.station +
        "</b><br>" +
        st.latitude.toFixed(4) +
        "\u00b0, " +
        st.longitude.toFixed(4) +
        "\u00b0<br>Elev: " +
        st.elevation.toFixed(0) +
        " m</div>",
      { className: "station-tooltip", direction: "top", offset: [0, -6] }
    );

    marker.bindPopup(
      '<div class="popup-title">' +
        st.network +
        "." +
        st.station +
        "</div>" +
        '<div class="popup-row"><span>Lat:</span> ' +
        st.latitude.toFixed(4) +
        "\u00b0</div>" +
        '<div class="popup-row"><span>Lon:</span> ' +
        st.longitude.toFixed(4) +
        "\u00b0</div>" +
        '<div class="popup-row"><span>Elev:</span> ' +
        st.elevation.toFixed(0) +
        " m</div>" +
        '<div class="popup-row"><span>Network:</span> ' +
        st.network +
        "</div>"
    );

    stationLayerGroup.addLayer(marker);
  });
}

// ---------------------------------------------------------------------------
// Event rendering
// ---------------------------------------------------------------------------
function magColor(mag) {
  if (mag === null || mag === undefined) return "#8890a8";
  if (mag < 2) return "#3ddc84";
  if (mag < 4) return "#ffb833";
  if (mag < 6) return "#ff8833";
  return "#f05050";
}

function magRadius(mag) {
  if (mag === null || mag === undefined) return 6;
  return Math.max(4, Math.min(25, 3 + mag * 3));
}

function magClass(mag) {
  if (mag === null || mag === undefined) return "low";
  if (mag < 3) return "low";
  if (mag < 5) return "mid";
  return "high";
}

function epochToUTC(epoch) {
  if (!epoch) return "—";
  const d = new Date(epoch * 1000);
  return d.toISOString().replace("T", " ").substring(0, 19) + " UTC";
}

function renderEvents(data) {
  events = data;
  eventLayerGroup.clearLayers();
  arrivalLayerGroup.clearLayers();

  const listEl = document.getElementById("event-list");
  listEl.innerHTML = "";

  if (!events || events.length === 0) {
    listEl.innerHTML = '<div style="padding:12px;color:var(--text-dim);font-size:12px;">No events detected yet.</div>';
    return;
  }

  events.forEach(function (ev) {
    // Map marker
    const color = magColor(ev.magnitude);
    const radius = magRadius(ev.magnitude);

    // Outer pulse ring for recent events
    const ageHours = (Date.now() / 1000 - ev.time) / 3600;
    if (ageHours < 1) {
      const pulse = L.circleMarker([ev.latitude, ev.longitude], {
        radius: radius,
        fillColor: color,
        fillOpacity: 0.2,
        color: color,
        weight: 1,
        opacity: 0.3,
        className: "event-pulse",
      });
      eventLayerGroup.addLayer(pulse);
    }

    const marker = L.circleMarker([ev.latitude, ev.longitude], {
      radius: radius,
      fillColor: color,
      fillOpacity: 0.7,
      color: "#fff",
      weight: 1.5,
    });

    const magStr =
      ev.magnitude !== null && ev.magnitude !== undefined
        ? (ev.magnitude_type || "M") + " " + ev.magnitude.toFixed(1)
        : "M ?";

    marker.bindPopup(
      '<div class="popup-title">Event ' +
        ev.evid +
        "</div>" +
        '<div class="popup-row"><span>Magnitude:</span> ' +
        magStr +
        "</div>" +
        '<div class="popup-row"><span>Lat:</span> ' +
        ev.latitude.toFixed(4) +
        "\u00b0</div>" +
        '<div class="popup-row"><span>Lon:</span> ' +
        ev.longitude.toFixed(4) +
        "\u00b0</div>" +
        '<div class="popup-row"><span>Depth:</span> ' +
        ev.depth.toFixed(1) +
        " km</div>" +
        '<div class="popup-row"><span>Time:</span> ' +
        epochToUTC(ev.time) +
        "</div>" +
        '<div class="popup-row"><span>Phases:</span> ' +
        (ev.num_phases || 0) +
        "</div>"
    );

    marker.on("click", function () {
      selectEvent(ev.evid);
    });

    eventLayerGroup.addLayer(marker);

    // Draw arrival lines
    if (ev.arrivals && ev.arrivals.length > 0) {
      ev.arrivals.forEach(function (arr) {
        const sta = stations.find(
          (s) => s.station === arr.station
        );
        if (sta) {
          const line = L.polyline(
            [
              [ev.latitude, ev.longitude],
              [sta.latitude, sta.longitude],
            ],
            {
              color: arr.phase === "P" ? "#4e9af1" : "#f0a050",
              weight: 1,
              opacity: 0.4,
              dashArray: "4 4",
            }
          );
          arrivalLayerGroup.addLayer(line);
        }
      });
    }

    // Sidebar list item
    const item = document.createElement("div");
    item.className = "event-item" + (ev.evid === selectedEventId ? " active" : "");
    item.dataset.evid = ev.evid;
    item.innerHTML =
      '<div class="ev-header">' +
      '<span class="ev-mag ' +
      magClass(ev.magnitude) +
      '">' +
      (ev.magnitude !== null && ev.magnitude !== undefined
        ? ev.magnitude.toFixed(1)
        : "?") +
      "</span>" +
      '<span class="ev-time">' +
      epochToUTC(ev.time).substring(0, 16) +
      "</span>" +
      "</div>" +
      '<div class="ev-loc">' +
      ev.latitude.toFixed(2) +
      "\u00b0, " +
      ev.longitude.toFixed(2) +
      "\u00b0 &middot; " +
      ev.depth.toFixed(1) +
      " km</div>" +
      '<div class="ev-phases">' +
      (ev.num_phases || 0) +
      " phases &middot; " +
      (ev.num_stations || 0) +
      " stations</div>";
    item.addEventListener("click", function () {
      selectEvent(ev.evid);
    });
    listEl.appendChild(item);
  });
}

function selectEvent(evid) {
  selectedEventId = evid;
  document.querySelectorAll(".event-item").forEach(function (el) {
    el.classList.toggle("active", el.dataset.evid == evid);
  });

  const ev = events.find((e) => e.evid == evid);
  if (ev) {
    map.flyTo([ev.latitude, ev.longitude], Math.max(map.getZoom(), 6), {
      duration: 0.8,
    });
  }
}

// ---------------------------------------------------------------------------
// Detail overlay
// ---------------------------------------------------------------------------
function showEventDetail(ev) {
  const magStr =
    ev.magnitude !== null && ev.magnitude !== undefined
      ? (ev.magnitude_type || "M") + " " + ev.magnitude.toFixed(2)
      : "—";

  document.getElementById("detail-title").textContent =
    "Event " + ev.evid;
  const tbody = document.getElementById("detail-body");
  tbody.innerHTML =
    "<tr><td>Origin Time</td><td>" +
    epochToUTC(ev.time) +
    "</td></tr>" +
    "<tr><td>Latitude</td><td>" +
    ev.latitude.toFixed(4) +
    "\u00b0</td></tr>" +
    "<tr><td>Longitude</td><td>" +
    ev.longitude.toFixed(4) +
    "\u00b0</td></tr>" +
    "<tr><td>Depth</td><td>" +
    ev.depth.toFixed(1) +
    " km</td></tr>" +
    "<tr><td>Magnitude</td><td>" +
    magStr +
    "</td></tr>" +
    "<tr><td>Phases</td><td>" +
    (ev.num_phases || 0) +
    "</td></tr>" +
    "<tr><td>Stations</td><td>" +
    (ev.num_stations || 0) +
    "</td></tr>";

  const arrDiv = document.getElementById("detail-arrivals");
  if (ev.arrivals && ev.arrivals.length > 0) {
    let html =
      '<table class="arrivals-table"><thead><tr>' +
      "<th>Station</th><th>Phase</th><th>Time</th><th>Residual</th><th>Dist</th><th>Az</th>" +
      "</tr></thead><tbody>";
    ev.arrivals.forEach(function (a) {
      html +=
        "<tr>" +
        "<td>" +
        a.station +
        "</td>" +
        "<td>" +
        (a.phase || "?") +
        "</td>" +
        "<td>" +
        (a.time ? a.time.toFixed(2) : "—") +
        "</td>" +
        "<td>" +
        (a.residual !== null && a.residual !== undefined ? a.residual.toFixed(3) : "—") +
        "</td>" +
        "<td>" +
        (a.distance !== null && a.distance !== undefined ? a.distance.toFixed(1) : "—") +
        "</td>" +
        "<td>" +
        (a.azimuth !== null && a.azimuth !== undefined ? a.azimuth.toFixed(1) : "—") +
        "\u00b0</td>" +
        "</tr>";
    });
    html += "</tbody></table>";
    arrDiv.innerHTML = html;
  } else {
    arrDiv.innerHTML =
      '<div style="color:var(--text-dim);font-size:12px;">No arrival data.</div>';
  }

  document.getElementById("event-detail").classList.remove("hidden");
}

function closeDetail() {
  document.getElementById("event-detail").classList.add("hidden");
}

// Allow double-click on event list items to open detail
document.getElementById("event-list").addEventListener("dblclick", function (e) {
  const item = e.target.closest(".event-item");
  if (item) {
    const ev = events.find((ev) => ev.evid == item.dataset.evid);
    if (ev) showEventDetail(ev);
  }
});

// ---------------------------------------------------------------------------
// Stats
// ---------------------------------------------------------------------------
function updateStats(stats) {
  document.getElementById("stat-stations").textContent =
    stats.station_count || 0;
  document.getElementById("stat-events").textContent =
    stats.event_count || 0;
  document.getElementById("stat-arrivals").textContent =
    stats.arrival_count || 0;
  document.getElementById("stat-origins").textContent =
    stats.origin_count || 0;
}

// ---------------------------------------------------------------------------
// Clock
// ---------------------------------------------------------------------------
function updateClock() {
  const now = new Date();
  document.getElementById("clock").textContent =
    now.toISOString().replace("T", " ").substring(0, 19) + " UTC";
}
setInterval(updateClock, 1000);
updateClock();

// ---------------------------------------------------------------------------
// WebSocket
// ---------------------------------------------------------------------------
function connectWebSocket() {
  const proto = location.protocol === "https:" ? "wss:" : "ws:";
  const url = proto + "//" + location.host + "/ws";

  ws = new WebSocket(url);

  ws.onopen = function () {
    document.getElementById("conn-status").textContent = "Connected";
    document.getElementById("conn-status").className =
      "status-badge connected";
    if (wsReconnectTimer) {
      clearTimeout(wsReconnectTimer);
      wsReconnectTimer = null;
    }
  };

  ws.onclose = function () {
    document.getElementById("conn-status").textContent = "Disconnected";
    document.getElementById("conn-status").className =
      "status-badge disconnected";
    scheduleReconnect();
  };

  ws.onerror = function () {
    document.getElementById("conn-status").textContent = "Error";
    document.getElementById("conn-status").className =
      "status-badge disconnected";
  };

  ws.onmessage = function (evt) {
    try {
      const msg = JSON.parse(evt.data);
      handleMessage(msg);
    } catch (e) {
      console.error("WS parse error:", e);
    }
  };
}

function scheduleReconnect() {
  if (!wsReconnectTimer) {
    wsReconnectTimer = setTimeout(function () {
      wsReconnectTimer = null;
      connectWebSocket();
    }, 3000);
  }
}

function handleMessage(msg) {
  switch (msg.type) {
    case "init":
      renderStations(msg.stations || []);
      renderEvents(msg.events || []);
      updateStats(msg.stats || {});
      fitMapToData();
      break;

    case "events_update":
      renderEvents(msg.events || []);
      if (msg.stats) updateStats(msg.stats);
      break;

    case "heartbeat":
      if (msg.stats) updateStats(msg.stats);
      break;

    case "pong":
      break;
  }
}

function fitMapToData() {
  const points = [];
  stations.forEach(function (s) {
    points.push([s.latitude, s.longitude]);
  });
  events.forEach(function (e) {
    points.push([e.latitude, e.longitude]);
  });
  if (points.length > 0) {
    const bounds = L.latLngBounds(points);
    map.fitBounds(bounds.pad(0.1));
  }
}

// Keep-alive ping
setInterval(function () {
  if (ws && ws.readyState === WebSocket.OPEN) {
    ws.send(JSON.stringify({ type: "ping" }));
  }
}, 30000);

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
document.addEventListener("DOMContentLoaded", function () {
  initMap();
  connectWebSocket();
});
